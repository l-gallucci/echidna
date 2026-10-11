"""Tests for the chained workflow commands (misassembly-wf, mag-quality-wf)
and their shared plan helpers. All assertions are on the pure build_plan() /
filter_plan() output, so no external tool is executed."""

import argparse
from pathlib import Path

import pytest

from mamisa.utils.workflow_run import Step, filter_plan, wrap_env
from mamisa.commands import misassembly_wf, mag_quality_wf, run_anvi_misassemblies


# ---- shared helpers -------------------------------------------------------

def _steps():
    return [Step(k, k, None, ['x']) for k in ['a', 'b', 'c', 'd']]


def test_wrap_env():
    assert wrap_env(None, ['mamisa', 'x']) == ['mamisa', 'x']
    assert wrap_env('gunc', ['mamisa', 'x']) == \
        ['conda', 'run', '--no-capture-output', '-n', 'gunc', 'mamisa', 'x']


def test_filter_plan_from_and_skip():
    s = _steps()
    assert [x.key for x in filter_plan(s, None, None)] == ['a', 'b', 'c', 'd']
    assert [x.key for x in filter_plan(s, 'c', None)] == ['c', 'd']
    assert [x.key for x in filter_plan(s, None, ['b'])] == ['a', 'c', 'd']
    assert [x.key for x in filter_plan(s, 'b', ['c'])] == ['b', 'd']


def test_filter_plan_unknown_key_returns_none():
    s = _steps()
    assert filter_plan(s, 'zzz', None) is None
    assert filter_plan(s, None, ['nope']) is None


# ---- misassembly-wf -------------------------------------------------------

def _misasm_args(**over):
    base = dict(output_dir=Path('OUT'), anvi_prefix='MisAsm',
                kraken2_output=None, kaiju_output=None,
                bam=Path('m.bam'), assembly=Path('asm.fa'), threads=8,
                anvio_env='anvio-9', checkm2_env='checkm2',
                max_length=300000, min_completeness=50.0, max_contamination=10.0,
                kaiju_names=None, kraken2_report=None)
    base.update(over)
    return argparse.Namespace(**base)


def test_misassembly_plan_without_taxonomy_skips_readchim():
    steps = misassembly_wf.build_plan(_misasm_args())
    keys = [s.key for s in steps]
    assert keys == ['anvi', 'process', 'classify', 'filter']
    # no chimera-report without read taxonomy
    filt = next(s for s in steps if s.key == 'filter')
    assert '--chimera-report' not in filt.cmd
    classify = next(s for s in steps if s.key == 'classify')
    assert '--taxonomy-windows' not in classify.cmd


def test_misassembly_plan_with_kraken_adds_readchim():
    steps = misassembly_wf.build_plan(
        _misasm_args(kraken2_output=Path('r.k2'), kraken2_report=Path('rep.txt')))
    keys = [s.key for s in steps]
    assert keys == ['anvi', 'process', 'readchim', 'classify', 'filter']
    rc = next(s for s in steps if s.key == 'readchim')
    assert '--kraken2-output' in rc.cmd and '--kraken2-report' in rc.cmd
    filt = next(s for s in steps if s.key == 'filter')
    assert '--chimera-report' in filt.cmd


def test_misassembly_anvi_runs_in_anvio_env_via_wrapper():
    steps = misassembly_wf.build_plan(_misasm_args())
    anvi = steps[0]
    assert anvi.env == 'anvio-9'
    assert anvi.cmd[:2] == ['mamisa', 'run-anvi-misassemblies']


# ---- mag-quality-wf -------------------------------------------------------

def _mag_args(**over):
    base = dict(output_dir=Path('MOUT'), genome_ext='fa', skip_mimag=False,
                bins_dir=Path('bins'), threads=40, checkm2_db=None,
                checkm2_env='checkm2', mamisa_env=None, tiers='HQ,MQ',
                gtdbtk_data=None, gtdbtk_env='gtdbtk', gunc_db=None,
                gunc_env='gunc', split_by_domain=False, mimag_env='mimag')
    base.update(over)
    return argparse.Namespace(**base)


def test_mag_plan_full_order():
    steps = mag_quality_wf.build_plan(_mag_args())
    assert [s.key for s in steps] == \
        ['checkm2', 'organize1', 'gtdbtk', 'gunc', 'mimag', 'chimeras', 'organize2']


def test_mag_plan_gunc_has_contig_taxonomy():
    steps = mag_quality_wf.build_plan(_mag_args(gunc_db=Path('db.dmnd')))
    gunc = next(s for s in steps if s.key == 'gunc')
    assert '--contig-taxonomy-output' in gunc.cmd
    assert '--db-file' in gunc.cmd


def test_mag_plan_gtdbtk_sets_data_path_env():
    steps = mag_quality_wf.build_plan(_mag_args(gtdbtk_data=Path('/data/gtdb')))
    gt = next(s for s in steps if s.key == 'gtdbtk')
    assert gt.env_vars == {'GTDBTK_DATA_PATH': '/data/gtdb'}


def test_mag_plan_skip_mimag_drops_step_and_gate():
    steps = mag_quality_wf.build_plan(_mag_args(skip_mimag=True))
    keys = [s.key for s in steps]
    assert 'mimag' not in keys
    org2 = next(s for s in steps if s.key == 'organize2')
    assert '--mimag-rna-dir' not in org2.cmd


def test_mag_plan_split_by_domain_wires_gtdbtk_dir():
    steps = mag_quality_wf.build_plan(_mag_args(split_by_domain=True))
    mi = next(s for s in steps if s.key == 'mimag')
    assert '--split-by-domain' in mi.cmd and '--gtdbtk-dir' in mi.cmd


def test_mag_check_chimeras_frozen_names_use_selected_dir():
    steps = mag_quality_wf.build_plan(_mag_args())
    chim = next(s for s in steps if s.key == 'chimeras')
    # check-chimeras reads the same Selected dir gunc saw (no rename between them)
    assert 'MOUT/02_filtered/Selected' in chim.cmd
    gunc = next(s for s in steps if s.key == 'gunc')
    assert 'MOUT/02_filtered/Selected' in gunc.cmd


# ---- run-anvi-misassemblies ----------------------------------------------

def test_anvi_just_do_it_default_on():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers()
    run_anvi_misassemblies.register_parser(sub)
    args = parser.parse_args(['run-anvi-misassemblies', '-b', 'x.bam',
                              '-o', 'pre/MisAsm'])
    assert args.just_do_it is True
    args2 = parser.parse_args(['run-anvi-misassemblies', '-b', 'x.bam',
                               '-o', 'pre/MisAsm', '--no-just-do-it'])
    assert args2.just_do_it is False
