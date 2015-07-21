#!/bin/bash
kill `ps -fu $USER|grep "master.R"|grep -v "grep"|cut -c10-14`
qdel -u $USER
