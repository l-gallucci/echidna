"""Shared helpers for the chained workflow commands (misassembly-wf, mag-quality-wf).

Each workflow is a list of Steps. A Step runs one command, optionally inside a
named conda env (`conda run -n <env> ...`) and optionally with extra environment
variables (e.g. GTDBTK_DATA_PATH). build_plan() in each workflow module is a pure
function so the exact command sequence can be unit-tested with --dry-run and
without any external tool installed.
"""

import os
import shlex
import subprocess
from typing import List, Optional, NamedTuple, Dict

from .logging import log_info, log_error, log_warning, print_section


class Step(NamedTuple):
    key: str                      # short id, usable with --from / --skip
    name: str                     # human label
    env: Optional[str]            # conda env to wrap in, or None for current env
    cmd: List[str]                # argv
    env_vars: Optional[Dict[str, str]] = None   # extra environment variables


def wrap_env(env: Optional[str], cmd: List[str]) -> List[str]:
    """Wrap a command in `conda run -n <env>` when env is set.

    --no-capture-output streams the tool's own stdout/stderr live, which matters
    for long-running steps (CheckM2, GTDB-Tk, GUNC, DIAMOND).
    """
    if env:
        return ['conda', 'run', '--no-capture-output', '-n', env] + cmd
    return cmd


def filter_plan(steps: List[Step], from_key: Optional[str],
                skip_keys: Optional[List[str]]) -> Optional[List[Step]]:
    """Apply --from (start at a step) and --skip (drop steps) to a plan.

    Returns None if from_key or any skip key is not a known step key (caller
    reports the error with the valid keys).
    """
    keys = [s.key for s in steps]
    skip = set(skip_keys or [])
    for k in list(skip) + ([from_key] if from_key else []):
        if k not in keys:
            return None
    start = keys.index(from_key) if from_key else 0
    return [s for i, s in enumerate(steps) if i >= start and s.key not in skip]


def run_step(idx: int, total: int, step: Step, dry_run: bool) -> bool:
    """Run one step. Returns True on success (or in dry-run)."""
    full = wrap_env(step.env, step.cmd)
    label = f"[{idx}/{total}] {step.name}"
    if step.env:
        label += f"  (env: {step.env})"
    print_section(label)
    prefix = ""
    if step.env_vars:
        prefix = " ".join(f"{k}={shlex.quote(v)}" for k, v in step.env_vars.items()) + " "
    log_info("  " + prefix + " ".join(shlex.quote(c) for c in full))

    if dry_run:
        return True

    run_env = None
    if step.env_vars:
        run_env = os.environ.copy()
        run_env.update(step.env_vars)

    try:
        subprocess.run(full, check=True, stdin=subprocess.DEVNULL, env=run_env)
        return True
    except subprocess.CalledProcessError as e:
        log_error(f"Step '{step.key}' failed (exit {e.returncode}). "
                  f"Fix it, then resume with --from {step.key}")
        return False
    except FileNotFoundError as e:
        log_error(f"Step '{step.key}': command not found ({e}). "
                  f"Is conda on PATH and the env created?")
        return False
    except KeyboardInterrupt:
        log_warning("\nInterrupted by user")
        return False


def execute_plan(steps: List[Step], from_key: Optional[str],
                 skip_keys: Optional[List[str]], dry_run: bool) -> int:
    """Filter then run a plan. Returns a process exit code (0 on success)."""
    plan = filter_plan(steps, from_key, skip_keys)
    if plan is None:
        valid = ', '.join(s.key for s in steps)
        log_error(f"Unknown step in --from/--skip. Valid steps: {valid}")
        return 2
    if not plan:
        log_warning("Nothing to run after applying --from/--skip")
        return 0

    total = len(plan)
    for i, step in enumerate(plan, 1):
        if not run_step(i, total, step, dry_run):
            return 1
    return 0
