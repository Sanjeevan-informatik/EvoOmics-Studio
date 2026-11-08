import gzip
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from evoomics.models import ValidationError
from evoomics.pipeline import run, sha256

ROOT = Path(__file__).resolve().parents[1]


class PipelineTests(unittest.TestCase):
    def test_end_to_end_and_hashes(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = run(ROOT / 'examples/aligned.fasta', Path(tmp) / 'run', model='jc69')
            manifest = json.loads((out / 'provenance.json').read_text())
            self.assertEqual(manifest['input']['sha256'], sha256(ROOT / 'examples/aligned.fasta'))
            for name, digest in manifest['artifacts_sha256'].items():
                self.assertEqual(sha256(out / name), digest)
            self.assertEqual(len(manifest['artifacts_sha256']), 7)
            self.assertTrue((out / 'tree.nwk').read_text().endswith(';\n'))
            with self.assertRaises(ValidationError):
                run(ROOT / 'examples/aligned.fasta', out)

    def test_repeatable_numerical_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            a = run(ROOT / 'examples/tiny.vcf', Path(tmp)/'a', input_format='vcf')
            b = run(ROOT / 'examples/tiny.vcf', Path(tmp)/'b', input_format='vcf')
            for name in ('tree.nwk', 'distances.tsv', 'qc.json'):
                self.assertEqual((a/name).read_bytes(), (b/name).read_bytes())

    def test_gzip(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / 'input.fasta.gz'
            with gzip.open(source, 'wb') as stream:
                stream.write((ROOT / 'examples/aligned.fasta').read_bytes())
            self.assertTrue(run(source, Path(tmp)/'out').exists())

    def test_snp_correction_guard(self):
        with tempfile.TemporaryDirectory() as tmp, self.assertRaisesRegex(ValidationError, 'ascertainment'):
            run(ROOT / 'examples/tiny.vcf', Path(tmp)/'out', input_format='vcf', model='jc69')

    def test_failure_leaves_no_run_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, out = Path(tmp)/'bad.fa', Path(tmp)/'out'
            source.write_text('>a\nNN\n>b\nNN\n')
            with self.assertRaises(ValidationError):
                run(source, out)
            self.assertFalse(out.exists())

    def test_cli_success_and_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            args = [sys.executable, '-m', 'evoomics', 'analyze', '--input',
                    str(ROOT/'examples/aligned.fasta'), '--format', 'fasta', '--output', str(Path(tmp)/'run')]
            result = subprocess.run(args, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            result = subprocess.run(args, capture_output=True, text=True)
            self.assertEqual(result.returncode, 2)
            self.assertNotIn('Traceback', result.stderr)
