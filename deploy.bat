@ECHO OFF
REM creates the zip file mosaic.zip with all files
REM needed for an installation of the application

ECHO running...

zip -r -9 mosaic.zip mods/*.py
zip -r -9 mosaic.zip data/*.ffa
zip -r -9 mosaic.zip doc/pic/*.gif
zip -r -9 mosaic.zip doc/pic/*.png
zip -r -9 mosaic.zip doc/papers
zip -r -9 mosaic.zip doc/index.html doc/style.css doc/manual.pdf
zip mosaic.zip example.bat mosaic.py readme.txt

ECHO finished.