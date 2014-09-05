#!/bin/bash
###############################################################################
#
#	Filename:	INSTALL_MIREVO (INSTALL)
#
#	Programmer: Ming Wen
#
#	$Id: INSTALL_MIREVO ,v 1.0 2009/10/06 
#
###############################################################################
#
#	Description:
#
#		This script installs the miREvo package
#
###############################################################################

# Make sure rigdir is defined appropriately
if [ $# -eq 0 ];then
	echo ""
	echo "    Usage: $0 dir_to_install_miREvo"
	echo ""
	exit 192;
fi

echo ""
echo ""

echo "Installing MIREVO ...... " 
export MIREVO=$1

echo "MIREVO is set to be \"$MIREVO\""

echo "Install miREvo to \"$MIREVO\""

echo ""
echo ""


if [ ! -e $MIREVO ]; then
	mkdir $MIREVO
fi
	
# Make sure we have the stuff we need
if [ -f miREvo/miREvo ]
then
	if [ -d $MIREVO ]
	then	
		if [ -w $MIREVO ]
		then
			# We have the proper permissions for $RIGDIR
			cp -r miREvo/* $MIREVO
			if [ $? -ne 0 ]
			then
				# We didn't succeed extracting the file!
				echo "Error in fetching miREvo"
				exit 1
			fi
		else
			echo "$MIREVO must be writable by the user installing the package."
			exit 1
		fi
	else
		echo "$MIREVO must be owned by the user installing the package."
		exit 1
	fi
else
	echo "Installation failed: "
	echo "A folder \"miREvo\" must exist in the current directory"
fi

echo "Installation succesed: "
echo ""
echo "Pleae set MIREVO to \"$MIREVO\" as an environment variable by export it in you .bashrc file, and also add it into you PERL runtime environment,"
echo ""
echo "Mostly this is done by adding: "
echo ""
echo "        export MIREVO=$MIREVO"
echo "        export PERL5LIB=\"\$PERL5LIB:\$MIREVO\""
echo ""
echo "into you ~/.bashrc file"
echo ""
echo "It is a good idea to put \"$MIREVO\" into your PATH after installing miREvo, by add the following line into you .barhrc file"
echo ""
echo "        export PATH=\"\$PATH:\$MIREVO\""
echo ""
echo "And use the command:" 
echo ""
echo "        source ~/.barhrc"
echo ""
echo "to make the above chagnes work."
echo ""
