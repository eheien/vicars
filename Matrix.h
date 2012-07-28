/*
 Top class representing a dense matrix for computation.
 */
class VCDenseMatrix {
protected:
  unsigned int  width, height;
public:
  VCDenseMatrix(const unsigned int &ncols, const unsigned int &nrows) : width(ncols), height(nrows) {};
  virtual void allocateRow(const unsigned int &row) = 0;
  virtual double val(const unsigned int &row, const unsigned int &col) const = 0;
  virtual void setVal(const unsigned int &row, const unsigned int &col, const double &new_val) = 0;
  virtual bool transpose(void) const = 0;
  virtual bool compressed(void) const = 0;
  virtual bool compressRow(const unsigned int &row, const float &ratio) = 0;
  virtual bool decompressRow(const unsigned int &row) = 0;
  virtual double *getRow(double *buf, const unsigned int &row) const = 0;
  virtual double *getCol(double *buf, const unsigned int &col) const = 0;
  virtual unsigned long mem_bytes(void) const = 0;
};

/*
 Subclass representing a non-compressed dense matrix.
 */
class VCDenseStd : public VCDenseMatrix {
protected:
  double   *data;
public:
  VCDenseStd(const int &ncols, const int &nrows);
  ~VCDenseStd(void);
  void allocateRow(const unsigned int &row) {};
  bool compressed(void) const { return false; };
  bool compressRow(const unsigned int &row, const float &ratio) { return false; };
  bool decompressRow(const unsigned int &row) { return false; };
  unsigned long mem_bytes(void) const;
};

/*
 Sub-subclass representing a normally oriented non-compressed dense matrix.
 */
class VCDenseStdStraight : public VCDenseStd {
public:
  VCDenseStdStraight(const unsigned int &ncols, const unsigned int &nrows) : VCDenseStd(ncols, nrows) {};
  bool transpose(void) const { return false; };
  double val(const unsigned int &row, const unsigned int &col) const { return data[row*width+col]; };
  void setVal(const unsigned int &row, const unsigned int &col, const double &new_val) { data[row*width+col] = new_val; };
  double *getRow(double *buf, const unsigned int &row) const;
  double *getCol(double *buf, const unsigned int &col) const;
};
