#!/bin/sh

MACRO=$1
let SOURCE_NUMBER=$2+1

sed -i s/SOURCE_NUMBER/$SOURCE_NUMBER/g $MACRO

WORK_DIR=`pwd`

export VO_CMS_SW_DIR=/uscmst1/prod/sw/cms
source $VO_CMS_SW_DIR/cmsset_default.sh
export SCRAM_ARCH=slc5_amd64_gcc462
scramv1 project CMSSW CMSSW_5_3_8_patch3

cd CMSSW_5_3_8_patch3
mv $WORK_DIR/src.tgz .
rm -r src/
tar -xzf src.tgz
cd src/SusyAnalysis/SusyNtuplizer/ttggAnalysis
eval `scramv1 runtime -sh`
make

mv $WORK_DIR/Cert_*_JSON*.txt .
mv $WORK_DIR/*.C .

root -b -q -l $MACRO

iscopy=1
numfail=0
while [[ $iscopy -ne 0 && $numfail -lt 10 ]]
do
	cp Run2012D-22Jan2013-DoublePhoton_$SOURCE_NUMBER\.root /eos/uscms/store/user/bfrancis/skims/ReReco/
        iscopy=$?
        numfail=$(( $numfail + 1 ))
done

echo
echo "Tried to copy $numfail times, ended with status $iscopy"
echo

if [ $iscopy -ne 0 ];
then
	mv Run2012D-22Jan2013-DoublePhoton_$SOURCE_NUMBER\.root $WORK_DIR
fi

cd $WORK_DIR
rm -rf CMSSW_5_3_8_patch3/
