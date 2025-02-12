# -*- coding: utf-8 -*-
"""
 * @Date: 2025-01-11 11:57:00
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2025-02-12 11:36:43
 * @FilePath: /pymummer/tests/test_genomediff.py
 * @Description:

 test:
    PYTHONPATH=. python tests/test_genomediff.py
 or:
    python -m pytest
"""
# """

from io import StringIO
from unittest import TestCase, main

from genomediff import Metadata, GenomeDiff
from genomediff.parser import GenomeDiffParser
from genomediff.records import Record


class ParserTestCase(TestCase):
    def test_parse(self):
        file = StringIO(
            "\n".join(
                [
                    "#=GENOME_DIFF\t1.0",
                    "#=AUTHOR test",
                    "SNP\t1\t23423\tNC_000913\t223\tA\tgene_name=mhpE",
                    "RA\t2\t\tNC_000913\t223\t0\tG\tA\tfrequency=0.1366",
                ]
            )
        )
        p = GenomeDiffParser(fsock=file)
        # fmt: off

        self.assertEqual(
            [
                Metadata("GENOME_DIFF", "1.0"),
                Metadata("AUTHOR", "test"),
                Record("SNP", 1, parent_ids=[23423], new_seq="A", seq_id="NC_000913", position=223, gene_name="mhpE"),
                Record("RA", 2, new_base="A", frequency=0.1366, position=223, seq_id="NC_000913", insert_position=0, ref_base="G"),
            ],
            list(p),
        )
        # fmt: on

    def test_parse_dot_missing_parent_ids(self):
        file = StringIO(
            "\n".join(
                [
                    "#=GENOME_DIFF\t1.0",
                    "#=AUTHOR test",
                    "SNP\t1\t23423\tNC_000913\t223\tA\tgene_name=mhpE",
                    "RA\t2\t.\tNC_000913\t223\t0\tG\tA\tfrequency=0.1366",
                ]
            )
        )
        p = GenomeDiffParser(fsock=file)
        # fmt: off
        self.assertEqual(
            [
                Metadata('GENOME_DIFF', '1.0'),
                Metadata('AUTHOR', 'test'),
                Record('SNP', 1, parent_ids=[23423], new_seq='A', seq_id='NC_000913', position=223, gene_name='mhpE'),
                Record('RA', 2, new_base='A', frequency=0.1366, position=223, seq_id='NC_000913', insert_position=0, ref_base='G')
            ],
            list(p)
        )
        # fmt: on


class GenomeDiffTestCase(TestCase):
    def test_document(self):
        file = StringIO(
            "\n".join(
                [
                    "#=GENOME_DIFF\t1.0",
                    "#=AUTHOR test",
                    "SNP\t1\t23423\tNC_000913\t223\tA",
                    "RA\t2\t\tNC_000913\t223\t0\tG\tA",
                ]
            )
        )
        document = GenomeDiff.read(file)

        # fmt: off
        self.assertEqual({'AUTHOR': ['test'], 'GENOME_DIFF': ['1.0']}, document.metadata._dict)
        self.assertEqual({'AUTHOR': 'test', 'GENOME_DIFF': '1.0'}, document.metadata)

        snp_record = Record('SNP', 1, [23423], seq_id='NC_000913', new_seq='A', position=223)
        ra_record = Record('RA', 2, None, position=223, seq_id='NC_000913', insert_position=0, new_base='A', ref_base='G')
        # fmt: on

        self.assertEqual([snp_record], document.mutations)
        self.assertEqual([ra_record], document.evidence)
        self.assertEqual(snp_record, document[1])
        self.assertEqual(ra_record, document[2])

    def test_query(self):
        file = StringIO(
            "\n".join(
                [
                    "#=GENOME_DIFF\t1.0",
                    "#=AUTHOR test",
                    "SNP\t1\t23423\tNC_000913\t223\tA",
                    "RA\t2\t\tNC_000913\t223\t0\tG\tA",
                ]
            )
        )
        document = GenomeDiff.read(file)

        snp_record = Record(
            "SNP", 1, [23423], seq_id="NC_000913", new_seq="A", position=223
        )
        ra_record = Record(
            "RA",
            2,
            None,
            position=223,
            seq_id="NC_000913",
            insert_position=0,
            new_base="A",
            ref_base="G",
        )

        self.assertEqual([snp_record], list(document.records.query(type="SNP")))
        self.assertEqual([snp_record], list(document.records.query(mut_type="SNP")))
        self.assertEqual([snp_record], list(document.records.query(type="==SNP")))
        self.assertEqual([snp_record], list(document.records.query("type==SNP")))
        self.assertEqual(
            [ra_record], list(document.records.query("position>0", mut_type="RA"))
        )


class RecordTestCase(TestCase):
    def test_simple(self):
        # fmt: off
        snp_record = Record('SNP', 1, parent_ids=[23423], seq_id='NC_000913', new_seq='A', position=223, test='more')

        self.assertEqual('SNP', snp_record.type)
        self.assertEqual(1, snp_record.id)
        self.assertEqual('A', snp_record.new_seq)
        self.assertEqual('more', snp_record.test)
        # fmt: on


class ParentResolveTestCase(TestCase):
    def test_resolve(self):
        file = StringIO(
            "\n".join(
                [
                    "#=GENOME_DIFF\t1.0",
                    "#=AUTHOR test",
                    "SNP\t1\t2\tNC_000913\t223\tA\tgene_name=mhpE",
                    "RA\t2\t\tNC_000913\t223\t0\tG\tA\tfrequency=0.1366",
                ]
            )
        )
        document = GenomeDiff.read(file)
        self.assertEqual(document.records.parents_of(1), [document[2]])
        self.assertEqual(document[1].parents, [document[2]])


class RecordComparisonTestCase(TestCase):
    def test_cmp1(self):
        file1 = StringIO(
            "\n".join(
                [
                    "#=GENOME_DIFF\t1.0",
                    "#=CREATED\t20:02:17 23 Jan 2019",
                    "#=PROGRAM\tbreseq 0.33.2",
                    "#=COMMAND\tbreseq -r LCA.gff3 sequence-data/DM0 evolved re-runs (Rohan)/ZDBp889_R1.fastq.gz sequence-data/DM0 evolved re-runs (Rohan)/ZDBp889_R2.fastq.gz sequence-data/ZDBp889_reads.fastq -o consensus/ZDBp889",
                    "#=REFSEQ\tLCA.gff3",
                    "#=READSEQ\tsequence-data/DM0 evolved re-runs (Rohan)/ZDBp889_R1.fastq.gz",
                    "#=READSEQ\tsequence-data/DM0 evolved re-runs (Rohan)/ZDBp889_R2.fastq.gz",
                    "#=READSEQ\tsequence-data/ZDBp889_reads.fastq",
                    "#=CONVERTED-BASES\t644779377",
                    "#=CONVERTED-READS\t14448149",
                    "#=INPUT-BASES\t645034321",
                    "#=INPUT-READS\t14455411",
                    "#=MAPPED-BASES\t602854657",
                    "#=MAPPED-READS\t13788351",
                    "SNP\t1\t34\tREL606\t72313\tC",
                ]
            )
        )
        document1 = GenomeDiff.read(file1)
        file2 = StringIO(
            "\n".join(
                [
                    "#=GENOME_DIFF\t1.0",
                    "#=CREATED\t16:49:49 23 Jan 2019",
                    "#=PROGRAM\tbreseq 0.33.2",
                    "#=COMMAND\tbreseq -r LCA.gff3 sequence-data/DM0 evolved re-runs (Rohan)/ZDB67_R1.fastq.gz sequence-data/DM0 evolved re-runs (Rohan)/ZDB67_R2.fastq.gz -o consensus/ZDB67",
                    "#=REFSEQ\tLCA.gff3",
                    "#=READSEQ\tsequence-data/DM0 evolved re-runs (Rohan)/ZDB67_R1.fastq.gz",
                    "#=READSEQ\tsequence-data/DM0 evolved re-runs (Rohan)/ZDB67_R2.fastq.gz",
                    "#=CONVERTED-BASES\t114566968",
                    "#=CONVERTED-READS\t419781",
                    "#=INPUT-BASES\t114567554",
                    "#=INPUT-READS\t419783",
                    "#=MAPPED-BASES\t92472620",
                    "#=MAPPED-READS\t339813",
                    "SNP\t1\t12\tREL606\t72313\tC",
                ]
            )
        )
        document2 = GenomeDiff.read(file2)
        self.assertEqual(document1.mutations, document2.mutations)

    def test_cmp2(self):
        file1 = StringIO(
            "\n".join(
                [
                    "#=GENOME_DIFF\t1.0",
                    "SNP\t1\t12\tREL606\t72313\tC\taa_new_seq=G\taa_position=92\taa_ref_seq=D\tcodon_new_seq=GGC\tcodon_number=92\tcodon_position=2\tcodon_ref_seq=GAC\tgene_name=araA\tgene_position=275\tgene_product=L-arabinose isomerase\tgene_strand=<\tgenes_overlapping=araA\tlocus_tag=ECB_00064\tlocus_tags_overlapping=ECB_00064\tmutation_category=snp_nonsynonymous\tposition_end=72313\tposition_start=72313\tsnp_type=nonsynonymous\ttransl_table=11",
                ]
            )
        )
        document1 = GenomeDiff.read(file1)
        file2 = StringIO(
            "\n".join(
                [
                    "#=GENOME_DIFF\t1.0",
                    "SNP\t1\t34\tREL606\t72313\tC\taa_new_seq=G\taa_position=92\taa_ref_seq=D\tcodon_new_seq=GGC\tcodon_number=92\tcodon_position=2\tcodon_ref_seq=GAC\tgene_name=araA\tgene_position=275\tgene_product=L-arabinose isomerase\tgene_strand=<\tgenes_overlapping=araA\tlocus_tag=ECB_00064\tlocus_tags_overlapping=ECB_00064\tmutation_category=snp_nonsynonymous\tposition_end=72313\tposition_start=72313\tsnp_type=nonsynonymous\ttransl_table=11",
                ]
            )
        )
        document2 = GenomeDiff.read(file2)
        self.assertEqual(document1.mutations, document2.mutations)


if __name__ == "__main__":
    main()
