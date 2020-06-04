# pdf2png.sh
# useful bit of code using imagemagick to convert pdf to png image
# from http://blog.robfelty.com/2008/03/11/convert-pdf-to-png-with-imagemagick/

for file in fig*; do \
echo $file;\
convert -density 600x600 -resize 1000x750 -quality 90 $file `echo $file|cut -f1 -d'.'`.png;\
done
