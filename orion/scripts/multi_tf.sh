#!/bin/sh
cd ~/Downloads/panelimp; W=orion/work
while [ ! -f $W/models.done ]; do sleep 30; done
/usr/bin/time -l target/release/panelimp multi --frame target --refs K562=$W/m_k562,HCT116=$W/m_hct116,HEK293T=$W/m_hek293t --targets RPE1=$W/m_rpe1,HepG2=$W/m_hepg2,Jurkat=$W/m_jurkat,K562=$W/m_k562,HCT116=$W/m_hct116,HEK293T=$W/m_hek293t --out $W/multi_tf > $W/multi_tf.log 2>&1
echo done > $W/multi_tf.done
