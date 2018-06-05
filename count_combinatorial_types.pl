# Extracts the different vertex-facets incidence graphs of K3 polytopes 
# with the same number of vertices.
# This is done by giving in input reflexives polytopes with the same number 
# of lattice points and their triangulations.

# The vertex-facets incidence graphs are stored in a text file in the form
# of the vertex-facets incidence matrix. 

# This has been used to extract the incidence graphs of K3 polytopes obtained
# from the reflexive polytopes with at most 30 lattice points in  
# "data/reflexive.txt" and their triangulations. 


use strict;
use warnings;
use application 'polytope';

# Sobstitute the first entry of each row such that each row sums up to 4
sub dehomogenize {
  my $m = new Matrix<Rational>($_[0]);
  for my $i (0..$m->rows-1){
    $m->elem($i,0) = 4 - $m->elem($i,1) - $m->elem($i,2) - $m->elem($i,3);
  }
  return $m;
}
 

# Given a matrix M with k rows, returns a map f:{0,...,k-1}->Rows(M)
sub point_function {
        my $matrix = $_[0];
        my $f = new Map<Int,Vector<Rational>>();
        for my $i (0..$matrix->rows-1) {
                $f->{$i} = $matrix->row($i);
        }
        return $f;
}

# Given a polytope P and a list of polytopes, checks whether the vertex-facets incidence graph is not isomorphic 
# to the ones of the polytopes contained in the list 
sub isnotcontained{
    my($poly, $list) = @_;
    if ($list->size==0){
	return true;
    }
    if ($list->size!=0){
	for my $elem (@{$list}){
	    if (isomorphic($poly->GRAPH->ADJACENCY,new Polytope(POINTS=>$elem)->GRAPH->ADJACENCY)==1){
		return false;
	    }
	}
    return true;
    }
}	   
		


my $set = new Set<Matrix>();
my $combinatorialtype = new Set<IncidenceMatrix<NonSymmetric>>();

open (my $file_all_reflexive, "<", "path_to_input_file_polytope.txt");
open (my $output, ">", "path_to_output_file");

my $i = 1;
my $k3 = new Polytope();
while(<$file_all_reflexive>){
    my $string = $_;
    my $M = new Matrix(eval $string);
    my $npts = $M->rows;
        open (my $file_triangulations, "<", "path_to_input_file_triangulation/".$i.".txt");
        while (<$file_triangulations>){
	    my $coefficients = $_;
	    my $hypersurface = new tropical::Hypersurface<Min>(MONOMIALS=>$M,COEFFICIENTS=>$coefficients);
	    for (my $k=0; $k <= $hypersurface->REGIONS->rows-1; $k++){
		my $p = new Polytope(POINTS=>$hypersurface->VERTICES->minor($hypersurface->REGIONS->[$k],All));
	        if ($p->BOUNDED==1){
		    $k3 = new Polytope(POINTS=>$p->VERTICES->minor(All,[0,2,3,4]));
		}
	    }
	    if (isnotcontained($k3, $set)==1){
		$set += $k3->VERTICES;
		$combinatorialtype += $k3->VERTICES_IN_FACETS;
	    }
	}
        ++$i;
}

close $file_all_reflexive;
print $output $combinatorialtype."\n";
print $output "number of combinatorial types:";
print $output $combinatorialtype->size;
close $output;
		   
