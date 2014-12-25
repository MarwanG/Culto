#!/bin/sh

current_dir=${0%/*}

cd $current_dir

echo $current_dir

PD_INST="/Applications/Pd-0.46-2.app/Contents/Resources"
PD_PATCH="sweep.pd"
PD_AUDIO="-r 44100 -audioindev 2 -audiooutdev 1 -audiobuf 46 -channels 2"
PD_MIDI="-nomidi"
PD_OPTIONS="-noprefs"
PD_PATH="-path ./abs:${PD_INST}/extra/iemabs::${PD_INST}/extra/zexy:${PD_INST}/extra/iemmatrix:${PD_INST}/extra/list-abs"
PD_LIB="-lib iemlib1:iemlib2:zexy:iemgui:iem_tab:iem_spec2:iem_ambi:iem_bin_ambi:iem_delay:iem_roomsim:iem_dp:iem_matrix:iemmatrix:net:OSC"


echo starting pd ...
${PD_INST}/bin/pd ${PD_AUDIO} ${PD_OPTIONS} ${PD_PATH} ${PD_LIB}
