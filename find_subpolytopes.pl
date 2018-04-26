# Find all the subpolytopes of the fourth dilation of the standard simplex
# up to permutation of the coordinates. In total, they are 356 461.
# The output is saved in the folder 'data/subpolytopes', sorted by number of
# lattice points.
#
# Among them, the reflexive polytopes are extracted and saved in the text file
# 'data/reflexive.txt'. They are 15 139.
#
# Finally, the minimal polytopes are extracted and saved in the text file
# 'data/minimal_reflexive.txt'. They are 115.

use strict;
use warnings;
use application 'polytope';

###############################################################################

# Sobstitute the first entry of each row of a matrix with 1
sub homogenize{
	my $m = new Matrix<Rational>($_[0]);
	for my $i (0..$m->rows-1){
		$m->elem($i,0) = 1;
	}
	return $m;
}

# Sobstitute the first entry of each row such that each row sums up to 4
sub dehomogenize{
	my $m = new Matrix<Rational>($_[0]);
	for my $i (0..$m->rows-1){
		$m->elem($i,0) = 4 - $m->elem($i,1) - $m->elem($i,2) - $m->elem($i,3);
	}
	return $m;
}

# Returns the lattice points of a polytope as a set
sub lattice_points{
	my $m = new Matrix<Rational>($_[0]->LATTICE_POINTS);
	my $s = new Set<Vector<Rational>>();
	for my $i (0..$m->rows-1){
		$s += $m->row($i);
	}
	return $s;
}

# Given a vector and a permutation, permutes its entries accordingly.
sub perm_coord{
	$v = $_[0];
	$p = $_[1];
	$new_v= new Vector<Rational>(1,$v->[($p->[0])+1],$v->[($p->[1])+1],$v->[($p->[2])+1],$v->[($p->[3])+1]);
	return $new_v;
}

###############################################################################

$PP = new Polytope(POINTS=>[[1,4,0,0,0],[1,0,4,0,0],[1,0,0,4,0],[1,0,0,0,4]]);

$f = new Map<Integer,Set<Set<Vector<Rational>>>>();
$f->{35} = new Set<Set<Vector<Rational>>>(lattice_points($PP));
for (my $size = 34; $size > 4; $size -= 1 ){
	my $empty_set = new Set<Set<Vector<Rational>>>();
	$f->{$size} = $empty_set;
}

for (my $size = 35; $size > 5; $size -= 1 ){
	for my $s_p (@{($f->{$size})}){
		my $p = new Polytope<Rational>(POINTS=>$s_p);
		my $verts_p = $p->VERTICES;
		for my $i (0..$verts_p->rows-1){
			my $s_q = $s_p - new Set<Vector<Rational>>($verts_p->row($i));	
			my $q = new Polytope(POINTS=>$s_q);
			if ($q->DIM eq 3 and $q->N_INTERIOR_LATTICE_POINTS eq 1){
				my $size_q = $q->N_LATTICE_POINTS;
				my $is_already_in = 0;
				for my $perm (@{group::all_group_elements(group::symmetric_group(4)->PERMUTATION_ACTION)}){
					my $perm_s_q = new Set<Vector<Rational>>();
					for my $v (@{$s_q}){
						$perm_s_q += perm_coord($v,$perm);
					}
					if (contains($f->{$size_q},$perm_s_q)){
						$is_already_in += 1;
						last;
					}
				}
				if ($is_already_in eq 0){
					$f->{$size_q} += $s_q;
				}
			}
		}
	}
	for (my $i = 35; $i > 4; $i -= 1 ){
		print $f->{$i}->size." ";
	}
	print "\n";
}

# Numbers of subpolytopes by number of lattice points (from 35 and decreasing)
# 1 1 2 5 14 30 76 168 365 727 1391 2474 4206 6753 10301 14922 20485 26589 32529 37336 39951 39654 36235 30179 22677 15040 8575 4014 1409 316 36 



