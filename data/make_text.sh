#!/bin/bash
set -eu

printf '' > lines.txt
printf '' > random.txt
printf '' > sorted.txt

for len in $(seq 1 5 40); do
    for i in $(seq 1 5000); do
        openssl rand -base64 $len >> lines.txt
    done
done

for i in $(seq 1 3); do
    shuf lines.txt >> random.txt
done

sort random.txt > sorted.txt