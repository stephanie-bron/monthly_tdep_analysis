
1. Make your own copy of labcore/start_lab.sh in the main lab directory:

  $ cd lab
  $ cp labcore/start_lab.sh  start.sh


2. Modify the user-variables in your start.sh so they point to your directories


3. Execute start.sh

  $ ./start.sh


4. Start Root

  $ root


5. Load the project you want (dependent projects load automatically), e.g:

  root[] .x llh/loadlibs.C

   Note that from any directory you can always load a library with e.g.:

  root[] .x $LAB_MAIN_DIR/llh/loadlibs.C


Trouble-shooting:

1. If you continue to have trouble loading or building a project, try
   deleting the corresponding lib directory, e.g.:

   $ cd $LAB_LIB_DIR
   $ rm -rf  llh  [or the name of whichever project is causing problems]

   Caution!  You are using -rf option; think twice before hitting enter
 
