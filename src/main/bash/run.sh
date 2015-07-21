#!/bin/bash
nohup lib/master.R $@ >screen_pipeline.out 2>screen_pipeline.err &
