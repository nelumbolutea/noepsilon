###bulk node child###
open(TR,"$ARGV[0]");
while(<TR>)
 {chomp;
  system("perl z_all_childnodes1.pl $_ > $_.nodes"); 
 }
####################