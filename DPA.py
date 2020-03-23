# Copyright (C) 2018 Daniel Page <csdsp@bristol.ac.uk>
#
# Use of this source code is restricted per the CC BY-NC-ND license, a copy of 
# which can be found via http://creativecommons.org (and should be included as 
# LICENSE.txt within the associated archive or repository).

import numpy, struct, sys, multiprocessing
from scipy.stats import pearsonr
from math import sqrt

sbox= [  
  0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
  0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
  0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
  0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
  0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
  0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
  0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
  0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
  0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
  0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
  0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
  0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
  0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
  0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
  0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
  0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16 ];



## Load  a trace data set from an on-disk file.
## 
## \param[in] f the filename to load  trace data set from
## \return    t the number of traces
## \return    s the number of samples in each trace
## \return    M a t-by-16 matrix of AES-128  plaintexts
## \return    C a t-by-16 matrix of AES-128 ciphertexts
## \return    T a t-by-s  matrix of samples, i.e., the traces

def traces_ld( f ) :
  fd = open( f, "rb" )

  def rd( x ) :
    ( r, ) = struct.unpack( x, fd.read( struct.calcsize( x ) ) ) ; return r

  t = rd( '<I' )
  s = rd( '<I' )

  M = numpy.zeros( ( t, 16 ), dtype = numpy.uint8 )
  C = numpy.zeros( ( t, 16 ), dtype = numpy.uint8 )
  T = numpy.zeros( ( t,  s ), dtype = numpy.int16 )

  for i in range( t ) :
    for j in range( 16 ) :
      M[ i, j ] = rd( '<B' )

  for i in range( t ) :
    for j in range( 16 ) :
      C[ i, j ] = rd( '<B' )

  for i in range( t ) :
    for j in range( s  ) :
      T[ i, j ] = rd( '<h' )

  fd.close()

  return t, s, M, C, T

## Store a trace data set into an on-disk file.
## 
## \param[in] f the filename to store trace data set into
## \param[in] t the number of traces
## \param[in] s the number of samples in each trace
## \param[in] M a t-by-16 matrix of AES-128  plaintexts
## \param[in] C a t-by-16 matrix of AES-128 ciphertexts
## \param[in] T a t-by-s  matrix of samples, i.e., the traces

def traces_st( f, t, s, M, C, T ) :
  fd = open( f, "wb" )

  def wr( x, y ) :
    fd.write( struct.pack( x, y ) )

  wr( '<I', t   )
  wr( '<I', s   )

  for i in range( t ) :
    for j in range( 16 ) :
      wr( '<B', M[ i, j ] )

  for i in range( t ) :
    for j in range( 16 ) :
      wr( '<B', C[ i, j ] )

  for i in range( t ) :
    for j in range( s  ) :
      wr( '<h', T[ i, j ] )

  fd.close()


def firstSbox(m,k):
  return sbox[m^k]

def hdist(x,y):
  return hweight(x^y)

#Hamming weight
def hweight(x):
  return bin(x).count('1')

def correlation(X,Y):
  x_bar = X.mean()
  y_bar = Y.mean()
  top = sum([(X[i]-x_bar) * (Y[i] - y_bar) for i in range(X.shape[0])])
  bottom = sqrt(sum([(x - x_bar)**2 for x in X]) * sum([(y-y_bar)**2 for y in Y]) )
  return(top/bottom)

## Attack implementation, as invoked from main after checking command line
## arguments.
##
## \param[in] argc number of command line arguments
## \param[in] argv           command line arguments

def attack( argc, argv ):
  t,s,M,C,T = traces_ld("traces.dat")
  M = M[:t]
  T = T[:t]
  K = numpy.array(range(0xff+1))
  print("Number of traces: {} \n Number of samples: {} \n Dim Messages: {} \n Dim Ciphertexts: {} \n Dim Traces: {}\n, Key Hypothesis {}".format(t,s,M.shape,C.shape,T.shape,K.shape))
  #The plan is to attack the the result of the first SBOX
  #So we create a hypothetical power consumption based on hamming leakage
  H = []
  for m in M:
    h = []
    for byte in m:
      h.append([hweight(byte,firstSbox(byte,k) ) for k in K])
    H.append(numpy.array(h))
  H = numpy.array(H)
  print("Dim H: {}, H[0] {}".format(H.shape,H[0]) )


  samples = T.shape[1] // 10
  with open("key.txt",'w') as f:
    for byte in range(M.shape[1]):
      print("Finding key for byte {}".format(byte))
      Hyp = H[:,byte,:]
      R = numpy.zeros((K.shape[0],samples))
      for t_index in range(samples):
        trace = T[:,t_index] # tc elements
        if t_index %1000 == 0:
          print("Reached t = ",t_index)
          print(numpy.unravel_index(numpy.argmax(R),R.shape), numpy.max(R) )
        for k in K:
          keyData = Hyp[:,k] #t elements
          cor, _ = pearsonr(keyData,trace)
          R[k,t_index] = 0 if numpy.isnan(cor) else abs(cor)
      print ("Likely Key: {} , Correlation {}".format(numpy.unravel_index(numpy.argmax(R),R.shape), numpy.max(R) ))
      f.write(str(numpy.unravel_index(numpy.argmax(R),R.shape)) + str(numpy.max(R)) )








if ( __name__ == '__main__' ) :
  attack( len( sys.argv ), sys.argv )