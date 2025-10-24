###keep original nodeID in first column###

use Bio::TreeIO;

my $treeio = Bio::TreeIO->new(-format => 'newick', -file => $ARGV[0]);
if (my $tree = $treeio->next_tree) {
    $n = 0;
    for my $node ($tree->get_nodes) {
        if ($node->is_Leaf) {
            
        } else {
            $n++;
            my $node_id = $node->id || "node$n"; 
            print "$node_id\t"; 

            for my $child ($node->each_Descendent) {
                my $child_id = $child->id;
                if ($child->is_Leaf) {
                    print $child_id . "\t";
                } else {
                    my $descendant = $tree->find_node(-id => $child_id);
                    if ($descendant->is_Leaf) {
                        print "\t";
                    } else {
                        for my $tip ($descendant->get_all_Descendents) {
                            if ($tip->is_Leaf) {
                                print $tip->id . ",";
                            }
                        }
                        print "\t";
                    }
                }
            }
            print "\n";
        }
    }
}

####