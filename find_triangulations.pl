# Extract all the central regular unimodular triangulations of the polytopes
# in a text file. This is done using the topcom implementation in 
# polymake. The triangulations are calculated on the boundary and then
# induced on the polytope.
#
# The triangulations are saved in a file (one for each polytope in the
# input text file). A single triangulation is stored as a list of heights
# such that the lower faces lifted polytope project to the triangulation.
#
#
# They have been calculated from the file 'data/reflexive.txt' for all the
# reflexive polytopes up to (normalised) volume 30. These are saved in the
# (compressed) folder 'data/triangulations/'
#
# They have been also calculated from the file 'data/minimal.txt' for all the
# minimal polytopes. These are saved in the folder 'data/triangulations_min/'

use strict;
use warnings;
use application 'polytope';

###############################################################################

# Given a matrix M with k rows, returns a map f:{0,...,k-1}->Rows(M)
sub point_function {
	my $matrix = $_[0];
	my $f = new Map<Int,Vector<Rational>>();
	for my $i (0..$matrix->rows-1) {
		$f->{$i} = $matrix->row($i);
	}
	return $f;
}

# Given a matrix with k rows, returns a map f:Rows(m)->{0,...,k-1}
sub index_function {
	my $matrix = $_[0];
	my $f = new Map<Vector<Rational>,Int>();
	for my $i (0..$matrix->rows-1) {
		$f->{$matrix->row($i)} = $i;
	}
	return $f;
}

# Returns the set of facets of a polytope, given as matrices
sub find_facets {
	my $polytope = $_[0];
	my $verts = $polytope->VERTICES;
	my $map_vert = point_function($verts);
	my $ViF = $polytope->VERTICES_IN_FACETS;
	my $faces = new Set<Matrix<Rational>>();
	for my $i (0..$ViF->rows-1){
		my $ind_verts = $ViF->row($i);
		my $face = new Matrix<Rational>(0,$polytope->DIM+1);
		for my $index (@{$ind_verts}){
			my $vertex = $map_vert->{$index};
			$face = $face/$vertex;
		}
		$faces += $face;
	}
	return $faces;
}

# Sobstitute the first entry of each row of a matrix with 1
sub homogenize {
	my $m = new Matrix<Rational>($_[0]);
	for my $i (0..$m->rows-1){
		$m->elem($i,0) = 1;
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



# Given a point configuration M, a subset m, and a triangulation t of m,
# relabel t with the indexing given by M. Also add the last point listed
# in M in each simplex of the new triangulation
sub relabel{
	my $M = $_[0];	# Matrix<Rational>
	my $m = $_[1];	# Matrix<Rational>
	my $t = $_[2];	# Array<Set<Int>>
	my $index_given_point = index_function($M);
	my $point_given_index = point_function($m);
	my $new_triangulation = new Set<Set<Int>>();
	for my $simplex (@{$t}){
		my $new_simplex = new Set<Int>();
		for my $index (@{$simplex}){
			$new_simplex += $index_given_point->{$point_given_index->{$index}};
		}
		$new_triangulation += $new_simplex;
	}
	return $new_triangulation;	# Set<Set<Int>>
}

# Given a set of sets as input, returns their cartesian product
sub cartesian {
	my $set=$_[0];
	my $all_poss = new Set<Set<Set<Int>>>();
	my $first = $set->[0];
	if ($set->size > 1){
		$set -= $first;
		for my $s (@{cartesian($set)}){
			for my $elem (@{$first}) {
				my $s2 = new Set<Set<Int>>();
				for (@{$elem}){
					$s2 += $_;
				}
				for (@{$s}){
					$s2 += $_;
				}
				$all_poss += $s2;
			}
		}
		return $all_poss;
	} else {
		for my $elem (@{$first}) {
			$all_poss += $elem;
		}
		return $all_poss;
	}
}

###############################################################################

open (my $input, "<", "path_to_input_file.txt");
my $i = 0;
my $tot = 0;
while(<$input>) {
	my $m = new Matrix<Rational>(eval $_);
	$i += 1;
	my $P = new Polytope<Rational>(POINTS=>homogenize($m));
	my $M = $P->BOUNDARY_LATTICE_POINTS;
	my $O = $P->INTERIOR_LATTICE_POINTS->row(0);
	my $pts_P = new Matrix<Rational>($M / $O);
	my $index_given_point = index_function($pts_P);
	my $n_triangs = 1;
	my $triangs_pyrs_relabeled = new Set<Set<Set<Set<Int>>>>();
	for my $facet (@{find_facets($P)}){
		my $verts_pyr_F = $facet/(new Vector<Rational>($O));
		my $pyr_F = new Polytope(VERTICES=>$verts_pyr_F);
		my $pts_pyr_F = new Matrix<Rational>($pyr_F->LATTICE_POINTS);
		my $pc_pyr_F = new PointConfiguration<Rational>(POINTS=>$pts_pyr_F);
		my $triangs_pyr_F = topcom_fine_and_connected_triangulations($pc_pyr_F);
		$n_triangs *= $triangs_pyr_F->size;
		# relabel the triangulation
		my $point_given_index = point_function($pts_pyr_F);
		my $triangs_pyr_F_relabeled = new Set<Set<Set<Int>>>();
		for my $triang (@{$triangs_pyr_F}){
			$triangs_pyr_F_relabeled += relabel($pts_P,$pts_pyr_F,$triang);
		}
		$triangs_pyrs_relabeled += $triangs_pyr_F_relabeled;
	}
	
		open (my $output, ">", "path_to_output_file/".$i.".txt") or die "Can't open : $!";
		my $triangs_P = cartesian($triangs_pyrs_relabeled);
		my $triangs_P_fixed_type = new Set<Array<Set<Int>>>();
		for (@{$triangs_P}){
			$triangs_P_fixed_type += new Array<Set<Int>>($_);
		}
		$triangs_P = $triangs_P_fixed_type;
		my $all_central_regular_triangs_P = new Set<Array<Set<Int>>>();
		for my $triang (@{$triangs_P}){
			#my $S = new fan::PolyhedralComplex(POINTS=>$pts_P,INPUT_POLYTOPES=>$triang);
			#$S->VISUAL;
			if (not is_subdivision($pts_P,$triang)){
				print "ERROR!\n";
				print $i;
				die;
			}
			my $pair = is_regular($pts_P,$triang);
			if ($pair->[0]){
				$all_central_regular_triangs_P += $triang;
				print $output $pair->[1]."\n";
				$tot += 1;
			}
		}
		print $i." - ".$all_central_regular_triangs_P->size."/".$triangs_P_fixed_type->size."\n";
		close $output;
}
close $input;
print $tot;



























