"""Run with python -m unittest src.analysis.mtec_high.simple_clustering.test_matching."""
import unittest
import polars as pl

from src.analysis.mtec_high.simple_clustering.helpers.match_clusters.clean_orthologs import clean_orthologs
from src.analysis.mtec_high.simple_clustering.helpers.match_clusters.match_mice_clusters import match_mice_clusters
from src.analysis.mtec_high.simple_clustering.helpers.match_clusters.match_rat_clusters import match_rat_clusters, select_shared_orthologs


class MatchingTests(unittest.TestCase):
    def setUp(self):
        self.frame = pl.DataFrame(
            [('ref', 'a', '1', True, 1), ('ref', 'b', '2', True, 2),
             ('ref', 'c', '3', True, None), ('ref', 'd', '4', True, 3),
             ('other', 'a', '1', True, 7), ('other', 'b', '2', True, 7)],
            schema=['dataset_id', 'gene_symbol', 'ortholog_id', 'has_clean_ortholog', 'local_cluster_id'],
            orient='row',
        )

    def test_gene_assignments_survive_split_clusters(self):
        mice = match_mice_clusters(self.frame, 'ref')
        self.assertEqual(mice.filter(pl.col('dataset_id') == 'other')['matched_cluster_id'].to_list(), [1, 2])
        rat = pl.DataFrame(
            [('rat_mtec_high', 'a', '1', True, 9), ('rat_mtec_high', 'b', '2', True, 9),
             ('rat_mtec_high', 'c', '3', True, None), ('rat_mtec_high', 'bad', '1///2', True, 10),
             ('rat_mtec_high', 'bad2', '5', False, 10)],
            schema=self.frame.schema, orient='row',
        )
        result = match_rat_clusters(mice, rat, 'ref')
        shared = select_shared_orthologs(result, 'ref')
        self.assertEqual(shared.filter(pl.col('dataset_id') == 'rat_mtec_high')['matched_cluster_id'].to_list(), [1, 2])
        self.assertEqual(shared['ortholog_id'].unique().sort().to_list(), ['1', '2'])

    def test_invalid_ids_and_duplicate_clean_ids(self):
        frame = pl.DataFrame(
            [('ref', '1', True), ('ref', 'A///B', True), ('ref', None, True),
             ('ref', ' ', True), ('ref', '2', False)],
            schema=['dataset_id', 'ortholog_id', 'has_clean_ortholog'], orient='row',
        )
        self.assertEqual(clean_orthologs(frame)['ortholog_id'].to_list(), ['1'])
        with self.assertRaises(ValueError):
            clean_orthologs(pl.concat([frame.head(1), frame.head(1)]))

    def test_rat_singleton_can_support_reference_region(self):
        mice = match_mice_clusters(self.frame, 'ref')
        rat = self.frame.head(1).with_columns(
            pl.lit('rat_mtec_high').alias('dataset_id'),
            pl.lit(None, dtype=pl.Int64).alias('local_cluster_id'),
        )
        result = select_shared_orthologs(match_rat_clusters(mice, rat, 'ref'), 'ref')
        self.assertEqual(result['matched_cluster_id'].unique().to_list(), [1])


if __name__ == '__main__':
    unittest.main()
