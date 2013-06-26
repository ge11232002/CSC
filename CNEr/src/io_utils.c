/****************************************************************************
 *                           I/O low-level utils                            *
 *                           Author: Ge Tan                                 *
 ****************************************************************************/
#include "CNEr.h"
#include  <ctype.h>  //for the tolower()
#include <sys/types.h>
#include <sys/stat.h> /* for fstat() */

/* --- .Call ENTRY POINT ---
 * Returns an external pointer.
 * From R:
 *   x <- .Call("new_input_ExternalFilePtr", "path/to/some/file",
 *              PACKAGE="Biostrings")
 *   reg.finalizer(x,
 *       function(e) .Call("ExternalFilePtr_close", e, PACKAGE="Biostrings"),
 *       onexit=TRUE)
 */
// Code taken from Biostrings

/* Code taken from do_url() in R/src/main/connections.c */
static void get_file_ztype(const char *path, int *ztype, int *subtype)
{
  FILE *fp;
  char buf[7];
  int res;

  *ztype = -1;
  *subtype = 0;
  if ((fp = fopen(path, "rb")) == NULL)
    return;
  memset(buf, 0, 7);
  res = fread(buf, 5, 1, fp);
  fclose(fp);
  if (res != 1)
    return;
  if (buf[0] == '\x1f' && buf[1] == '\x8b')
    *ztype = 0;
  else if (strncmp(buf, "BZh", 3) == 0)
    *ztype = 1;
  else if (buf[0] == '\xFD' && strncmp(buf+1, "7zXZ", 4) == 0)
    *ztype = 2;
  else if ((buf[0] == '\xFF') && strncmp(buf+1, "LZMA", 4) == 0) {
    *ztype = 2;
    *subtype = 1;
  } else if (memcmp(buf, "]\0\0\200\0", 5) == 0) {
    *ztype = 2;
    *subtype = 1;
  }
  return;
}

static FILE *open_file(const char *path, const char *mode)
{
  FILE *fp;
  int ret;
  struct stat buf;

  fp = fopen(path, mode);
  if (fp == NULL)
    error("cannot open file '%s'", path);
  ret = fstat(fileno(fp), &buf);
  if (ret != 0) {
    fclose(fp);
    error("Biostrings internal error in open_input_file(): "
          "cannot stat file '%s'", path);
  }
  if (S_ISDIR(buf.st_mode)) {
    fclose(fp);
    error("file '%s' is a directory", path);
  }
  return fp;
}

static FILE *open_input_file(const char *path)
{
  FILE *fp;
  int ret, ztype, subtype;
  struct stat buf;

  get_file_ztype(path, &ztype, &subtype);
  switch (ztype) {
  /* No compression */
      case -1:
    fp = open_file(path, "r");
      break;
  /* gzfile */
      case 0:
    error("cannot open file '%s' (gzip-compressed files "
          "are not supported yet, sorry!)", path);
    //open the file with gzFile gfp = gzopen(path, "r");
    //close it with gzclose();
      break;
  /* bzfile */
      case 1:
    error("cannot open file '%s' (bzip2-compressed files "
          "are not supported yet, sorry!)", path);
    //requires #include <bzlib.h>
    //open the file with BZFILE* bfp = BZ2_bzReadOpen(...);
    //close it with BZ2_bzReadClose();
      break;
  /* xzfile */
      case 2:
    error("cannot open file '%s' (LZMA-compressed files "
          "are not supported yet, sorry!)", path);
    //requires #include <lzma.h>
    //opening/closing this type of file seems quite complicated
      break;
      default:
    error("Biostrings internal error in open_input_file(): ",
          "invalid ztype value %d", ztype);
  }
  return fp;
}

SEXP new_input_ExternalFilePtr(SEXP filepath)
{
  SEXP filepath_elt;
  const char *expath;
  FILE *fp;

  if (!IS_CHARACTER(filepath) || LENGTH(filepath) != 1)
    error("'filepath' must be a single string");
  filepath_elt = STRING_ELT(filepath, 0);
  if (filepath_elt == NA_STRING)
    error("'filepath' is NA");
  // Maybe do this in R, not here.
  expath = R_ExpandFileName(translateChar(filepath_elt));
  fp = open_input_file(expath);
  return R_MakeExternalPtr(fp, R_NilValue, R_NilValue);
}

/* --- .Call ENTRY POINT ---
 * Closes the file pointed to by e.
 */
SEXP ExternalFilePtr_close(SEXP e)
{
  if (R_ExternalPtrAddr(e) != NULL) {
    //Rprintf("closing file...\n");
    fclose(R_ExternalPtrAddr(e));
    R_SetExternalPtrAddr(e, NULL);
  }
  return R_NilValue;
}



