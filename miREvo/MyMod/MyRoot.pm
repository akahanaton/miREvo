#
#===============================================================================
#
#         FILE:  MyRoot.pm
#
#  DESCRIPTION:  A Root class of my modules
#
#        FILES:  MyMod/MyRoot.pm
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Shen Yang 
#      COMPANY:  
#      VERSION:  1.0
#      CREATED:  03/24/2006 03:22:13 AM CST
#     REVISION:  ---
#===============================================================================
package MyMod::MyRoot;
use vars qw(@ISA $Default_Source);
use strict;
use warnings;
use Bio::Root::Root;
@ISA=qw(Bio::Root::Root);

sub _load_error{
	my ($self,$module)=@_;
	print STDERR <<END;
$self: $module cannot be found
Exception $@
for more information about the SeqIO system please see the SeqIO docs.
This includes ways of checking for formats at compile time, not run time
END
;
return 1;
}

sub _init_param{
	my ($self,@args)=@_;
	my %param=@args;
  @param{ map { lc $_ } keys %param } = values %param; # lowercase keys
  return %param;
}
