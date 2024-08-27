#!/bin/bash

grep -c ">" * | awk '{print $2 " " $1}'
