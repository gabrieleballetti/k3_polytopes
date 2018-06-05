# Given a reflexive polytope and a central regular unimodular triangulation, computes the K3 polytope
# defined by the qurtic surface dual to the triangulation.

#This has been used to analyze the K3 polytopes obtained from the minimal polytopes "data/minimal" and their triangulations

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



# Given as input monomials and coefficients, computes the tropical surfaces and extract the bounded region

sub compute_K3polytope{
    my($j, $dehom_monomials) = @_;
    open (my $input, "<", "path_to_input_file_triangulations/".$j.".txt");
    while (<$input>) {
	 my $coefficients = $_;
	 my $hypersurface = new tropical::Hypersurface<Min>(MONOMIALS=>$dehom_monomials,COEFFICIENTS=>$coefficients);
	 for (my $k=0; $k <= $hypersurface->REGIONS->rows-1; $k++){
	     my $p = new Polytope(POINTS=>$hypersurface->VERTICES->minor($hypersurface->REGIONS->[$k],All));
	     if ($p->BOUNDED==1) {
	         my $q = new Polytope(POINTS=>$p->VERTICES->minor(All,[0,2,3,4]));
		 return $q;		 
	     }
	 }
    }
    close $input;
}


open (my $polytopes, "<", "path_to_input_file_polytopes.txt");
open (my $output, ">", "path_to_output_file.txt");

my $i = 1;
while(<$polytopes>){
    print "i: $i\n";
    my $line = $_;
    my $monomials = new Matrix(eval $line);
    my $dehom_monomials = dehomogenize($monomials);
    my $K3 = compute_K3polytope($i, $dehom_monomials);
    print $output  $K3->F_VECTOR"\n";
    print $output $K3->VERTICES_IN_FACETS"\n";
    # $K3->VISUAL;
    ++$i;
}

close $polytopes;
close $output;

