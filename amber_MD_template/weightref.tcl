mol load pdb ref_raw.pdb

set all [atomselect top "all"]
$all set occupancy 0
$all set beta 0

set prot [atomselect top "noh protein"]
set lig [atomselect top "noh resname LIG"]

$prot set beta 0
$lig set beta 1
$prot set occupancy 1
$lig set occupancy 0

$all writepdb ref.pdb
quit
