for og in test/*; do
  ./diff-output/diff-output $og ./$(basename $og)
done
