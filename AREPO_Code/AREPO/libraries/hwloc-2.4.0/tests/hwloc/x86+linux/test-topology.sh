#!/bin/bash
#-*-sh-*-

#
# Copyright © 2015-2020 Inria.  All rights reserved.
# See COPYING in top-level directory.
#

HWLOC_top_srcdir="/home/astro/Desktop/Astro-and-Cosmo-Simulations/AREPO_Code/AREPO/libraries/hwloc-2.4.0"
HWLOC_top_builddir="/home/astro/Desktop/Astro-and-Cosmo-Simulations/AREPO_Code/AREPO/libraries/hwloc-2.4.0"
lstopo="$HWLOC_top_builddir/utils/lstopo/lstopo-no-graphics"

HWLOC_PLUGINS_PATH=${HWLOC_top_builddir}/hwloc/.libs
export HWLOC_PLUGINS_PATH

HWLOC_DONT_ADD_VERSION_INFO=1
export HWLOC_DONT_ADD_VERSION_INFO

HWLOC_XML_EXPORT_SUPPORT=0
export HWLOC_XML_EXPORT_SUPPORT

actual_output="$1"

# make sure we use default numeric formats
LANG=C
LC_ALL=C
export LANG LC_ALL

error()
{
    echo $@ 2>&1
}

# test_topology NAME TOPOLOGY-DIR
#
# Test the topology under TOPOLOGY-DIR.  Return true on success.
test_topology ()
{
    local name="$1"
    local dir="$2"
    local expected_output="$3"
    local options="$4"

    local output="`mktemp`"

    export HWLOC_THISSYSTEM=0
    export HWLOC_DEBUG_CHECK=1
    export HWLOC_COMPONENTS=x86,linux,stop
    export HWLOC_FSROOT="$dir/fsroot"
    export HWLOC_CPUID_PATH="$dir/cpuid"

    opts="--of xml -"
    [ -r "$options" ] && opts=`cat $options`

    if ! "$lstopo" $opts \
	| sed	-e 's/ gp_index="[0-9]*"//' \
	> "$output"
	# filtered gp_index because it may change if we reorder discovery
	# (not used in NUMA/PU distances, and regenerated (differently) during XML import).
    then
	result=1
    else
	if [ "$HWLOC_UPDATE_TEST_TOPOLOGY_OUTPUT" != 1 ]
	then
	    /usr/bin/diff -b -u "$expected_output" "$output"
	    result=$?
	else
	    if ! /usr/bin/diff "$expected_output" "$output" >/dev/null
	    then
		cp -f "$output" "$expected_output"
		echo "Updated $expected_output"
	    fi
	    result=0
	fi
    fi

    rm "$output"

    return $result
}

# test_eligible TOPOLOGY-DIR
#
# Return true if the topology under TOPOLOGY-DIR is eligible for
# testing with the current flavor.
test_eligible()
{
    local dir="$1"
    local output="$2"

    [ -d "$dir" -a -f "$output" ]
}


if [ ! -x "$lstopo" ]
then
    error "Could not find executable file \`$lstopo'."
    exit 1
fi

topology="${actual_output%.output}"
if [ "$topology" = "$actual_output" ] ;
then
    error "Input file \`$1' should end with .output"
    exit 1
fi
actual_options="$topology".options

# if there's a .source file, use the tarball name it contains instead of $topology
if [ -f "$topology".source ] ; then
    actual_source="$HWLOC_top_srcdir"/tests/hwloc/x86+linux/`cat "$topology".source`
else
    actual_source="$topology".tar.bz2
fi

# if there's a .env file, source it
if [ -f "$topology".env ] ; then
    . "$topology".env
fi

result=1

dir="`mktemp -d`"

if ! ( bunzip2 -c "$actual_source" | ( cd "$dir" && tar xf - $tar_options ) )
then
    error "failed to extract x86 cpuid \`$topology'"
else
    actual_dir="`echo "$dir"/*`"

    if test_eligible "$actual_dir" "$actual_output"
    then
	test_count="`expr $test_count + 1`"

	test_topology "`basename $topology`" "$actual_dir" "$actual_output" "$actual_options"
	result=$?
    else
	# Skip this test.
	result=77
    fi
fi

rm -rf "$dir"

exit $result