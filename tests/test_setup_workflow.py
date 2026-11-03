"""Tests for setup-workflow / check-envs pure helpers."""

from echidna.commands.setup_workflow import set_yaml_scalar
from echidna.commands.check_envs import _parse_ver, check_env


def test_set_yaml_scalar_replaces_existing():
    text = "threads: 1\nfoo: bar\n"
    out = set_yaml_scalar(text, "threads", 40)
    assert "threads: 40" in out
    assert "foo: bar" in out            # other lines preserved


def test_set_yaml_scalar_quotes():
    out = set_yaml_scalar("genomes_dir: x\n", "genomes_dir", "/a b/c", quote=True)
    assert 'genomes_dir: "/a b/c"' in out


def test_set_yaml_scalar_appends_when_missing():
    out = set_yaml_scalar("threads: 1\n", "env_gunc", "gunc", quote=True)
    assert 'env_gunc: "gunc"' in out


def test_set_yaml_scalar_only_top_level():
    # an indented key of the same name must not be touched
    text = "block:\n  threads: 1\nthreads: 1\n"
    out = set_yaml_scalar(text, "threads", 8)
    assert "  threads: 1" in out        # nested untouched
    assert "threads: 8" in out


def test_parse_ver():
    assert _parse_ver("samtools 1.21") == (1, 21, 0)
    assert _parse_ver("diamond version 2.1.24") == (2, 1, 24)
    assert _parse_ver("no numbers here") is None


def test_check_env_unreachable(monkeypatch):
    # when conda run fails, env check returns not ok without raising
    import echidna.commands.check_envs as ce
    monkeypatch.setattr(ce, "_conda_run", lambda env, cmd, timeout=120: (1, "no env"))
    result = ce.check_env("ghost", "gunc", repo="/repo")
    assert result["ok"] is False
    assert result["env"] == "ghost"


def test_patch_env_distribution_and_checkout_with_spaces(tmp_path):
    import json
    import shlex
    from echidna.commands.setup_workflow import patch_env_files
    envs = tmp_path / 'envs'
    envs.mkdir()
    env = envs / 'test.yaml'
    env.write_text('dependencies:\n  - pip:\n      - echidnaxomics\n')
    repo = tmp_path / 'repo with spaces'
    assert patch_env_files(envs, repo, False) == 1
    entry = json.loads(env.read_text().splitlines()[-1].strip()[2:])
    assert shlex.split(entry) == ['-e', str(repo)]
    before = env.read_text()
    assert patch_env_files(envs, repo, False) == 1
    assert env.read_text() == before


def test_export_installed_workflow(tmp_path, monkeypatch):
    import argparse
    import pytest
    from echidna.commands import setup_workflow as sw
    package = tmp_path / 'site-packages' / 'echidna'
    templates = package / 'workflow_templates'
    (templates / 'envs').mkdir(parents=True)
    (templates / 'config.yaml').write_text('threads: 1\n')
    (templates / 'Snakefile').write_text('# fixture\n')
    (templates / 'envs' / 'cli.yaml').write_text('  - echidnaxomics\n')
    monkeypatch.setattr(sw, 'PACKAGE_ROOT', package)
    parser = argparse.ArgumentParser()
    sw.register_parser(parser.add_subparsers())
    dest = tmp_path / 'exported'
    args = parser.parse_args(['setup-workflow', '--export-workflow', str(dest), '--threads', '2'])
    sw.run(args)
    assert 'threads: 2' in (dest / 'config.yaml').read_text()
    assert 'echidnaxomics' in (dest / 'envs' / 'cli.yaml').read_text()
    with pytest.raises(SystemExit):
        sw.run(args)  # an existing directory is never overwritten by export


def test_export_dry_run_does_not_write(tmp_path, monkeypatch):
    import argparse
    from echidna.commands import setup_workflow as sw
    package = tmp_path / 'site-packages' / 'echidna'
    (package / 'workflow_templates').mkdir(parents=True)
    monkeypatch.setattr(sw, 'PACKAGE_ROOT', package)
    parser = argparse.ArgumentParser()
    sw.register_parser(parser.add_subparsers())
    dest = tmp_path / 'exported'
    sw.run(parser.parse_args(['setup-workflow', '--export-workflow', str(dest), '--dry-run']))
    assert not dest.exists()
