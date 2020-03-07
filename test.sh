for x in out/*; do
  echo $x
  diff-output/diff-output $x data/$(basename $x)
done
