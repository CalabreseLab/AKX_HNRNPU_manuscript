#!/usr/bin/env bash

shopt -s nullglob
for neg in *_neg_readcounts.txt; do
  base=${neg%_neg_readcounts.txt}            
  pos=${base}_pos_readcounts.txt             
  out=${base}_readcounts.txt                 

  # read the two numbers (handles integers or floats)
  a=$(<"$neg")
  b=$(<"$pos")

  # compute sum (use bc for floats; if all integers you can do: sum=$((a + b)))
  sum=$((a + b))

  # write it out
  echo "$sum" > "$out"
done
