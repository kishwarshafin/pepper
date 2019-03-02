from build import HELEN
from graphviz import Digraph

from modules.python.Options import DeBruijnGraphOptions


class DeBruijnHaplotyper:
    def __init__(self, fasta_handler, contig, start, end):
        self.fasta_handler = fasta_handler
        self.contig = contig
        self.region_start = start
        self.region_end = end

    @staticmethod
    def visualize(graph, output_filename, on_pruned_nodes=False):
        dot = Digraph(comment='dbruijn graph')
        for i in range(1, graph.current_hash_value):
            if on_pruned_nodes is True and i not in graph.good_nodes:
                continue
            kmer = graph.node_hash_int_to_str[i]
            dot.node(str(i), label=kmer)

        for i in range(1, graph.current_hash_value):
            if on_pruned_nodes is True and i not in graph.good_nodes:
                continue
            if i not in graph.out_nodes:
                continue
            for j in range(len(graph.out_nodes[i])):
                node_a = i
                node_b = graph.out_nodes[i][j]
                if on_pruned_nodes is True and node_b not in graph.good_nodes:
                    continue
                weight = graph.edges[(node_a, node_b)]
                dot.edge(str(node_a), str(node_b), label=str(weight))

        print("GRAPH SAVED")
        # print(self.dot.source)
        dot.render('outputs/'+output_filename+'.sv')

    def find_haplotypes(self, reads):
        # get the reference from the fasta file
        reference_sequence = self.fasta_handler.get_reference_sequence(self.contig, self.region_start, self.region_end)
        min_k, max_k = FRIDAY.DeBruijnGraph.find_min_k_from_ref(reference_sequence,
                                                                DeBruijnGraphOptions.MIN_K,
                                                                DeBruijnGraphOptions.MAX_K,
                                                                DeBruijnGraphOptions.STEP_K)
        # couldn't build ref without cycle
        if min_k == -1:
            return None, None

        # print(reference_sequence)
        # min_k = 44
        for kmer_size in range(min_k, max_k+1, DeBruijnGraphOptions.STEP_K):
            dbg_graph = FRIDAY.DeBruijnGraph(self.region_start, self.region_end)
            haplotypes = dbg_graph.generate_haplotypes(reference_sequence, reads, kmer_size)
            # self.visualize(dbg_graph, 'unpruned', False)
            # break
            if haplotypes:
                # print(kmer_size)
                # self.visualize(dbg_graph, 'pruned', True)
                return reference_sequence, haplotypes

        return reference_sequence, []
