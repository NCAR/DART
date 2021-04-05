<!--

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

-->
<html>
<head>
<META http-equiv="Content-Type" content="text/html; charset=EUC-JP">
<title>Module fft_mod</title>
<link type="text/css" href="http://www.gfdl.noaa.gov/~fms/style/doc.css" rel="stylesheet">
<STYLE TYPE="text/css">
          .fixed {
            font-size:medium;
            font-family:monospace;
            border-style:none;
            border-width:0.1em;
            padding:0.1em;
            color:#663366;
          }
        </STYLE>
</head>
<body>
<a name="TOP"></a><font class="header" size="1"><a href="#PUBLIC INTERFACE">PUBLIC INTERFACE </a>~
          <a href="#PUBLIC DATA">PUBLIC DATA </a>~
          <a href="#PUBLIC ROUTINES">PUBLIC ROUTINES </a>~
          <a href="#NAMELIST">NAMELIST </a>~
          <a href="#DIAGNOSTIC FIELDS">DIAGNOSTIC FIELDS </a>~
          <a href="#ERROR MESSAGES">ERROR MESSAGES </a>~
          <a href="#REFERENCES">REFERENCES </a>~ 
          <a href="#NOTES">NOTES</a></font>
<hr>
<h2>module fft_mod</h2>
<a name="HEADER"></a>
<!-- BEGIN HEADER -->
<div>
<b>Contact:</b>&nbsp;<a href="mailto:bw@gfdl.noaa.gov">   Bruce Wyman </a>
<br>
<b>Reviewers:</b>&nbsp;<br>
<b>Change History: </b>&nbsp;<a href="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/">WebCVS Log</a>
<br>
<b>Last Modified:</b>&nbsp;2002/03/22 00:10:54<br>
<br>
</div>
<!-- END HEADER -->
<a name="OVERVIEW"></a>
<hr>
<h4>OVERVIEW</h4>
<!-- BEGIN OVERVIEW -->
<p class="text">   Performs simultaneous fast Fourier transforms (FFTs) between
   real grid space and complex Fourier space. </p>
<!-- END OVERVIEW -->
<a name="DESCRIPTION"></a>
<!-- BEGIN DESCRIPTION -->
<div>   This routine computes multiple 1-dimensional FFTs and inverse FFTs.
   There are 2d and 3d versions between type real grid point space
   and type complex Fourier space. There are single (32-bit) and
   full (64-bit) versions.
   <br>
<br>
   On Cray and SGI systems, vendor-specific scientific library
   routines are used, otherwise a user may choose a NAG library version
   or stand-alone version using Temperton's FFT. </div>
<br>
<!-- END DESCRIPTION -->
<a name="OTHER MODULES USED"></a>
<hr>
<h4>OTHER MODULES USED</h4>
<!-- BEGIN OTHER MODULES USED -->
<div>
<pre>platform_mod<br>     fms_mod<br>   fft99_mod</pre>
</div>
<!-- END OTHER MODULES USED -->
<!-- BEGIN PUBLIC INTERFACE -->
<a name="PUBLIC INTERFACE"></a>
<hr>
<h4>PUBLIC INTERFACE</h4>
<div>
<pre>use fft_mod [, only:  fft_grid_to_fourier,<br>                      fft_fourier_to_grid,<br>                      fft_init,<br>                      fft_end ]</pre>
<dl>
<dt>
<a href="#fft_grid_to_fourier">fft_grid_to_fourier</a>:</dt>
<dd>   Given multiple sequences of real data values, this routine
   computes the complex Fourier transform for all sequences. </dd>
<dt>
<a href="#fft_fourier_to_grid">fft_fourier_to_grid</a>:</dt>
<dd>   Given multiple sequences of Fourier space transforms,
   this routine computes the inverse transform and returns
   the real data values for all sequences. </dd>
<dt>
<a href="#fft_init">fft_init</a>:</dt>
<dd>   This routine must be called to initialize the size of a
   single transform and setup trigonometric constants. </dd>
<dt>
<a href="#fft_end">fft_end</a>:</dt>
<dd>   This routine is called to unset the transform size and deallocate memory. </dd>
</dl>
</div>
<br>
<!-- END PUBLIC INTERFACE -->
<a name="PUBLIC DATA"></a>
<hr>
<h4>PUBLIC DATA</h4>
<!-- BEGIN PUBLIC DATA -->
<div>None.<br>
<br>
</div>
<!-- END PUBLIC DATA -->
<a name="PUBLIC ROUTINES"></a>
<hr>
<h4>PUBLIC ROUTINES</h4>
<!-- BEGIN PUBLIC ROUTINES -->
<ol type="a">
<li>
<a name="fft_grid_to_fourier"></a>
<h4>fft_grid_to_fourier</h4>
<pre>fourier = <b>fft_grid_to_fourier</b> ( grid )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Given multiple sequences of real data values, this routine
   computes the complex Fourier transform for all sequences. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>grid&nbsp;&nbsp;&nbsp;</tt></td><td>   Multiple sequence of real data values. The first dimension
   must be n+1 (where n is the size of a single sequence). <br>&nbsp;&nbsp;&nbsp;<span class="type">[real(R4_KIND), dimension(:,:)]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[real(R8_KIND), dimension(:,:)]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[real(R4_KIND), dimension(:,:,:)]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[real(R8_KIND), dimension(:,:,:)]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>fourier&nbsp;&nbsp;&nbsp;</tt></td><td>   Multiple sequences of transformed data in complex Fourier space.
   The first dimension must equal n/2+1 (where n is the size
   of a single sequence). The remaining dimensions must be the
   same size as the input argument "grid". <br>&nbsp;&nbsp;&nbsp;<span class="type">[complex(R4_KIND), dimension(lenc,size(grid,2))]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[complex(R8_KIND), dimension(lenc,size(grid,2))]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[complex(R4_KIND), dimension(lenc,size(grid,2),size(grid,3))]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[complex(R8_KIND), dimension(lenc,size(grid,2),size(grid,3))]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>NOTE</b>
</dt>
<dd>   The complex Fourier components are passed in the following format. <pre>        fourier (1)     = cmplx ( a(0), b(0) )
        fourier (2)     = cmplx ( a(1), b(1) )
            :              :
            :              :
        fourier (n/2+1) = cmplx ( a(n/2), b(n/2) )</pre>   where n = length of each real transform </dd>
<br>
<br>
</dl>
</li>
<li>
<a name="fft_fourier_to_grid"></a>
<h4>fft_fourier_to_grid</h4>
<pre>grid = <b>fft_fourier_to_grid</b> ( fourier )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Given multiple sequences of Fourier space transforms,
   this routine computes the inverse transform and returns
   the real data values for all sequences. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>fourier&nbsp;&nbsp;&nbsp;</tt></td><td>   Multiple sequence complex Fourier space transforms.
   The first dimension must equal n/2+1 (where n is the
   size of a single real data sequence). <br>&nbsp;&nbsp;&nbsp;<span class="type">[real(R4_KIND), dimension(:,:)]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[real(R8_KIND), dimension(:,:)]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[real(R4_KIND), dimension(:,:,:)]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[real(R8_KIND), dimension(:,:,:)]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>grid&nbsp;&nbsp;&nbsp;</tt></td><td>   Multiple sequence of real data values. The first dimension
   must be n+1 (where n is the size of a single sequence).
   The remaining dimensions must be the same size as the input
   argument "fourier". <br>&nbsp;&nbsp;&nbsp;<span class="type">[complex(R4_KIND), dimension(leng1,size(fourier,2))]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[complex(R8_KIND), dimension(leng1,size(fourier,2))]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[complex(R4_KIND), dimension(leng1,size(fourier,2),size(fourier,3))]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[complex(R8_KIND), dimension(leng1,size(fourier,2),size(fourier,3))]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="fft_init"></a>
<h4>fft_init</h4>
<pre>
<b>call fft_init </b>( n )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   This routine must be called once to initialize the size of a
   single transform. To change the size of the transform the
   routine fft_exit must be called before re-initialing with fft_init. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>n&nbsp;&nbsp;&nbsp;</tt></td><td>   The number of real values in a single sequence of data.
   The resulting transformed data will have n/2+1 pairs of
   complex values. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="fft_end"></a>
<h4>fft_end</h4>
<pre>
<b>call fft_end </b>
</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   This routine is called to unset the transform size and
   deallocate memory. It can not be called unless fft_init
   has already been called. There are no arguments. </dd>
<br>
<br>
</dl>
</li>
</ol>
<!-- END PUBLIC ROUTINES -->
<a name="NAMELIST"></a>
<!-- BEGIN NAMELIST -->
<!-- END NAMELIST --><a name="DIAGNOSTIC FIELDS"></a>
<!-- BEGIN DIAGNOSTIC FIELDS -->
<!-- END DIAGNOSTIC FIELDS --><a name="DATA SETS"></a>
<!-- BEGIN DATA SETS -->
<hr>
<h4>DATA SETS</h4>
<div>None.<br>
<br>
</div>
<!-- END DATA SETS -->
<a name="ERROR MESSAGES"></a>
<!-- BEGIN ERROR MESSAGES -->
<hr>
<h4>ERROR MESSAGES</h4>
<div>
<dl>
<dt>
<b>Error in fft_grid_to_fourier</b>
</dt>
<dd>
<span class="errmsg">fft_init must be called</span>
</dd>
<dd>   The initialization routine fft_init must be called before routines
   fft_grid_to_fourier. </dd>
<dt>
<b>Error in fft_grid_to_fourier</b>
</dt>
<dd>
<span class="errmsg">size of first dimension of input data is wrong</span>
</dd>
<dd>   The real grid point field must have a first dimension equal to n+1
   (where n is the size of each real transform). This message occurs
   when using the SGI/Cray fft. </dd>
<dt>
<b>Error in fft_grid_to_fourier</b>
</dt>
<dd>
<span class="errmsg">length of input data too small</span>
</dd>
<dd>   The real grid point field must have a first dimension equal to n
   (where n is the size of each real transform). This message occurs
   when using the NAG or Temperton fft. </dd>
<dt>
<b>Error in fft_grid_to_fourier</b>
</dt>
<dd>
<span class="errmsg">float kind not supported for nag fft</span>
</dd>
<dd>   32-bit real data is not supported when using the NAG fft. You
   may try modifying this part of the code by uncommenting the
   calls to the NAG library or less consider using the Temperton fft. </dd>
<dt>
<b>Error in fft_fourier_to_grid</b>
</dt>
<dd>
<span class="errmsg">fft_init must be called</span>
</dd>
<dd>   The initialization routine fft_init must be called before routines fft_fourier_to_grid. </dd>
<dt>
<b>Error in fft_fourier_to_grid</b>
</dt>
<dd>
<span class="errmsg">size of first dimension of input data is wrong</span>
</dd>
<dd>   The complex Fourier field must have a first dimension equal to
   n/2+1 (where n is the size of each real transform). This message
   occurs when using the SGI/Cray fft. </dd>
<dt>
<b>Error in fft_fourier_to_grid</b>
</dt>
<dd>
<span class="errmsg">length of input data too small</span>
</dd>
<dd>   The complex Fourier field must have a first dimension greater
   than or equal to n/2+1 (where n is the size of each real
   transform). This message occurs when using the NAG or Temperton fft. </dd>
<dt>
<b>Error in fft_fourier_to_grid</b>
</dt>
<dd>
<span class="errmsg">float kind not supported for nag fft</span>
</dd>
<dd>   float kind not supported for nag fft 
   32-bit real data is not supported when using the NAG fft. You
   may try modifying this part of the code by uncommenting the
   calls to the NAG library or less consider using the Temperton fft. </dd>
<dt>
<b>FATAL in fft_init</b>
</dt>
<dd>
<span class="errmsg">attempted to reinitialize fft</span>
</dd>
<dd>   You must call fft_exit before calling fft_init for a second time. </dd>
<dt>
<b>Error in fft_end</b>
</dt>
<dd>
<span class="errmsg">attempt to un-initialize fft that has not been initialized</span>
</dd>
<dd>   You can not call fft_end unless fft_init has been called. </dd>
</dl>
<br>
</div>
<!-- END ERROR MESSAGES -->
<a name="REFERENCES"></a>
<hr>
<h4>REFERENCES</h4>
<!-- BEGIN REFERENCES -->
<div>
<ol>
<li>   For the SGI/Cray version refer to the manual pages for
   DZFFTM, ZDFFTM, SCFFTM, and CSFFTM. </li>
<li>   For the NAG version refer to the NAG documentation for
   routines C06FPF, C06FQF, and C06GQF. </li>
</ol>
</div>
<br>
<!-- END REFERENCES -->
<a name="COMPILER SPECIFICS"></a>
<hr>
<h4>COMPILER SPECIFICS</h4>
<!-- BEGIN COMPILER SPECIFICS -->
<div>
        None.
      </div>
<br>
<!-- END COMPILER SPECIFICS -->
<a name="PRECOMPILER OPTIONS"></a>
<hr>
<h4>PRECOMPILER OPTIONS</h4>
<!-- BEGIN PRECOMPILER OPTIONS -->
<div>
<dl>
<dt>-D<b> NAGFFT</b>
</dt>
<dd>   -D NAGFFT
   On non-Cray/SGI machines, set to use the NAG library FFT routines.
   Otherwise the Temperton FFT is used by default. </dd>
<dt>-D<b> test_fft</b>
</dt>
<dd>   Provides source code for a simple test program.
   The program generates several sequences of real data.
   This data is transformed to Fourier space and back to real data,
   then compared to the original real data. </dd>
</dl>
</div>
<br>
<!-- END PRECOMPILER OPTIONS -->
<a name="LOADER OPTIONS"></a>
<hr>
<h4>LOADER OPTIONS</h4>
<!-- BEGIN LOADER -->
<div>
<p>   On SGI machines the scientific library needs to be loaded by
   linking with: </p>
<pre>        -lscs</pre>
<p>   If using the NAG library, the following loader options (or
   something similar) may be necessary: </p>
<pre>        -L/usr/local/lib -lnag</pre>
</div>
<!-- END LOADER OPTIONS -->
<a name="TEST PROGRAM"></a>
<hr>
<h4>TEST PROGRAM</h4>
<!-- BEGIN TEST PROGRAM -->
<div>None.<br>
</div>
<br>
<!-- END TEST PROGRAM -->
<a name="KNOWN BUGS"></a>
<hr>
<h4>KNOWN BUGS</h4>
<!-- BEGIN KNOWN BUGS -->
<div>
        None.
      </div>
<br>
<!-- END KNOWN BUGS -->
<a name="NOTES"></a>
<hr>
<h4>NOTES</h4>
<!-- BEGIN NOTES -->
<div>   The routines are overloaded for 2d and 3d versions.
   The 2d versions copy data into 3d arrays then calls the 3d interface.
   <br>
<br>
   On SGI/Cray machines:
   <br>
<br>
   There are single (32-bit) and full (64-bit) versions.
   For Cray machines the single precision version does not apply.
   <br>
<br>
   On non-SGI/CRAY machines:
   <br>
<br>
   The NAG library option uses the "full" precision NAG
   routines (C06FPF,C06FQF,C06GQF). Users may have to specify
   a 64-bit real compiler option (e.g., -r8).
   <br>
<br>
   The stand-alone Temperton FFT option works for the
   real precision specified at compile time.
   If you compiled with single (32-bit) real precision
   then FFT's cannot be computed at full (64-bit) precision. </div>
<br>
<!-- END NOTES -->
<a name="FUTURE PLANS"></a>
<hr>
<h4>FUTURE PLANS</h4>
<!-- BEGIN FUTURE PLANS -->
<div>
        None.
      </div>
<br>
<!-- END FUTURE PLANS -->
<hr>
<div align="right">
<font size="-1"><a href="#TOP">top</a></font>
</div>
</body>
</html>
