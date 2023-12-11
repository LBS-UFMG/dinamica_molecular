gmx trjconv -s md.tpr -f md.xtc -b 0 -e 10000 -o md_noPBC.xtc -pbc mol -center -ur compact << eof
1
0
eof
gmx rms -s md.tpr -f md_noPBC.xtc -b 0 -e 10000 -o rmsd.xvg -tu ns << eof
4
4
eof
echo 4 | gmx rmsf -s md.tpr -f md_noPBC.xtc -b 0 -e 10000 -o rmsf_residue.xvg -res
