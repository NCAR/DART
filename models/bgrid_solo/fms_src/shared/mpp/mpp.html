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
<title>Module mpp_mod</title>
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
<h2>module mpp_mod</h2>
<a name="HEADER"></a>
<!-- BEGIN HEADER -->
<div>
<b>Contact:&nbsp;</b><a href="mailto:vb@gfdl.noaa.gov">   V. Balaji </a>
<br>
<b>Reviewers:&nbsp;</b>
<br>
<b>Change History:&nbsp;</b><a href="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/">WebCVS Log</a>
<br>
<b>RCS Log:&nbsp;</b><a href="http://www.gfdl.noaa.gov/~vb/changes_models/bgrid_solo/fms_src/shared/mpp/mpp.html">RCS Log</a>
<br>
<b>Last Modified:&nbsp;</b>2002/07/19 18:10:05<br>
<br>
</div>
<!-- END HEADER -->
<a name="OVERVIEW"></a>
<hr>
<h4>OVERVIEW</h4>
<!-- BEGIN OVERVIEW -->
<p class="text"> 
<tt>mpp_mod</tt>, is a set of simple calls to provide a uniform interface
   to different message-passing libraries. It currently can be
   implemented either in the SGI/Cray native SHMEM library or in the MPI
   standard. Other libraries (e.g MPI-2, Co-Array Fortran) can be
   incorporated as the need arises. </p>
<!-- END OVERVIEW -->
<a name="DESCRIPTION"></a>
<!-- BEGIN DESCRIPTION -->
<div>   The data transfer between a processor and its own memory is based
   on <tt>load</tt>   and <tt>store</tt>   operations upon
   memory. Shared-memory systems (including distributed shared memory
   systems) have a single address space and any processor can acquire any
   data within the memory by <tt>load</tt>   and <tt>store</tt>. The situation is different for distributed
   parallel systems. Specialized MPP systems such as the T3E can simulate
   shared-memory by direct data acquisition from remote memory. But if
   the parallel code is distributed across a cluster, or across the Net,
   messages must be sent and received using the protocols for
   long-distance communication, such as TCP/IP. This requires a
   ``handshaking'' between nodes of the distributed system. One can think
   of the two different methods as involving <tt>put</tt>s or <tt>get</tt>s (e.g the SHMEM library), or in the case of
   negotiated communication (e.g MPI), <tt>send</tt>s and <tt>recv</tt>s.
   <br>
<br>
   The difference between SHMEM and MPI is that SHMEM uses one-sided
   communication, which can have very low-latency high-bandwidth
   implementations on tightly coupled systems. MPI is a standard
   developed for distributed computing across loosely-coupled systems,
   and therefore incurs a software penalty for negotiating the
   communication. It is however an open industry standard whereas SHMEM
   is a proprietary interface. Besides, the <tt>put</tt>s or <tt>get</tt>s on which it is based cannot currently be implemented in
   a cluster environment (there are recent announcements from Compaq that
   occasion hope).
   <br>
<br>
   The message-passing requirements of climate and weather codes can be
   reduced to a fairly simple minimal set, which is easily implemented in
   any message-passing API. <tt>mpp_mod</tt>   provides this API.
   <br>
<br>
   Features of <tt>mpp_mod</tt>   include:
   <br>
<br>
   1) Simple, minimal API, with free access to underlying API for
   more complicated stuff.<br>   2) Design toward typical use in climate/weather CFD codes.<br>   3) Performance to be not significantly lower than any native API.
   <br>
<br>
   This module is used to develop higher-level calls for <a href="models/bgrid_solo/fms_src/shared/mpp/mpp_domains.html">domain decomposition</a>   and <a href="models/bgrid_solo/fms_src/shared/mpp/mpp_io.html">parallel I/O</a>.
   <br>
<br>
   Parallel computing is initially daunting, but it soon becomes
   second nature, much the way many of us can now write vector code
   without much effort. The key insight required while reading and
   writing parallel code is in arriving at a mental grasp of several
   independent parallel execution streams through the same code (the SPMD
   model). Each variable you examine may have different values for each
   stream, the processor ID being an obvious example. Subroutines and
   function calls are particularly subtle, since it is not always obvious
   from looking at a call what synchronization between execution streams
   it implies. An example of erroneous code would be a global barrier
   call (see <a href="#mpp_sync"></a>   below) placed
   within a code block that not all PEs will execute, e.g:
   <br>
<br> 
<pre>   if( pe.EQ.0 )call mpp_sync()</pre>   Here only PE 0 reaches the barrier, where it will wait
   indefinitely. While this is a particularly egregious example to
   illustrate the coding flaw, more subtle versions of the same are
   among the most common errors in parallel code.
   <br>
<br>
   It is therefore important to be conscious of the context of a
   subroutine or function call, and the implied synchronization. There
   are certain calls here (e.g <tt>mpp_declare_pelist, mpp_init,
   mpp_malloc, mpp_set_stack_size</tt>) which must be called by all
   PEs. There are others which must be called by a subset of PEs (here
   called a <tt>pelist</tt>) which must be called by all the PEs in the <tt>pelist</tt>   (e.g <tt>mpp_max, mpp_sum, mpp_sync</tt>). Still
   others imply no synchronization at all. I will make every effort to
   highlight the context of each call in the MPP modules, so that the
   implicit synchronization is spelt out.  
   <br>
<br>
   For performance it is necessary to keep synchronization as limited
   as the algorithm being implemented will allow. For instance, a single
   message between two PEs should only imply synchronization across the
   PEs in question. A <i>global</i>   synchronization (or <i>barrier</i>)
   is likely to be slow, and is best avoided. But codes first
   parallelized on a Cray T3E tend to have many global syncs, as very
   fast barriers were implemented there in hardware.
   <br>
<br>
   Another reason to use pelists is to run a single program in MPMD
   mode, where different PE subsets work on different portions of the
   code. A typical example is to assign an ocean model and atmosphere
   model to different PE subsets, and couple them concurrently instead of
   running them serially. The MPP module provides the notion of a <i>current pelist</i>, which is set when a group of PEs branch off
   into a subset. Subsequent calls that omit the <tt>pelist</tt>   optional
   argument (seen below in many of the individual calls) assume that the
   implied synchronization is across the current pelist. The calls <tt>mpp_root_pe</tt>   and <tt>mpp_npes</tt>   also return the values
   appropriate to the current pelist. The <tt>mpp_set_current_pelist</tt>   call is provided to set the current pelist. </div>
<br>
<!-- END DESCRIPTION -->
<a name="OTHER MODULES USED"></a>
<hr>
<h4>OTHER MODULES USED</h4>
<!-- BEGIN OTHER MODULES USED -->
<div>
<pre>shmem_interface<br>            mpi</pre>
</div>
<!-- END OTHER MODULES USED -->
<!-- BEGIN PUBLIC INTERFACE -->
<a name="PUBLIC INTERFACE"></a>
<hr>
<h4>PUBLIC INTERFACE</h4>
<div>   F90 is a strictly-typed language, and the syntax pass of the
   compiler requires matching of type, kind and rank (TKR). Most calls
   listed here use a generic type, shown here as <tt>MPP_TYPE_</tt>. This
   is resolved in the pre-processor stage to any of a variety of
   types. In general the MPP operations work on 4-byte and 8-byte
   variants of <tt>integer, real, complex, logical</tt>   variables, of
   rank 0 to 5, leading to 48 specific module procedures under the same
   generic interface. Any of the variables below shown as <tt>MPP_TYPE_</tt>   is treated in this way. <pre>use mpp_mod [, only:  mpp_max,<br>                      mpp_sum,<br>                      mpp_transmit,<br>                      mpp_broadcast,<br>                      mpp_chksum,<br>                      mpp_error,<br>                      mpp_init,<br>                      stdin,<br>                      mpp_exit,<br>                      mpp_pe,<br>                      mpp_npes,<br>                      mpp_declare_pelist,<br>                      mpp_set_current_pelist,<br>                      mpp_clock_set_grain,<br>                      mpp_sync,<br>                      mpp_sync_self,<br>                      mpp_malloc,<br>                      mpp_set_stack_size ]</pre>
<dl>
<dt>
<a href="#mpp_max">mpp_max</a>:</dt>
<dd>   Reduction operations. </dd>
<dt>
<a href="#mpp_sum">mpp_sum</a>:</dt>
<dd>   Reduction operation. </dd>
<dt>
<a href="#mpp_transmit">mpp_transmit</a>:</dt>
<dd>   Basic message-passing call. </dd>
<dt>
<a href="#mpp_broadcast">mpp_broadcast</a>:</dt>
<dd>   Parallel broadcasts. </dd>
<dt>
<a href="#mpp_chksum">mpp_chksum</a>:</dt>
<dd>   Parallel checksums. </dd>
<dt>
<a href="#mpp_error">mpp_error</a>:</dt>
<dd>   Error handler. </dd>
<dt>
<a href="#mpp_init">mpp_init</a>:</dt>
<dd>   Initialize <tt>mpp_mod</tt>. </dd>
<dt>
<a href="#stdin">stdin</a>:</dt>
<dd>   Standard fortran unit numbers. </dd>
<dt>
<a href="#mpp_exit">mpp_exit</a>:</dt>
<dd>   Exit <tt>mpp_mod</tt>. </dd>
<dt>
<a href="#mpp_pe">mpp_pe</a>:</dt>
<dd>   Returns processor ID. </dd>
<dt>
<a href="#mpp_npes">mpp_npes</a>:</dt>
<dd>   Returns processor count for current pelist. </dd>
<dt>
<a href="#mpp_declare_pelist">mpp_declare_pelist</a>:</dt>
<dd>   Declare a pelist. </dd>
<dt>
<a href="#mpp_set_current_pelist">mpp_set_current_pelist</a>:</dt>
<dd>   Set context pelist. </dd>
<dt>
<a href="#mpp_clock_set_grain">mpp_clock_set_grain</a>:</dt>
<dd>   Set the level of granularity of timing measurements. </dd>
<dt>
<a href="#mpp_sync">mpp_sync</a>:</dt>
<dd>   Global synchronization. </dd>
<dt>
<a href="#mpp_sync_self">mpp_sync_self</a>:</dt>
<dd>   Local synchronization. </dd>
<dt>
<a href="#mpp_malloc">mpp_malloc</a>:</dt>
<dd>   Symmetric memory allocation. </dd>
<dt>
<a href="#mpp_set_stack_size">mpp_set_stack_size</a>:</dt>
<dd>   Allocate module internal workspace. </dd>
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
<a name="mpp_max"></a>
<h4>mpp_max</h4>
<pre>
<b>call mpp_max </b>( a, pelist )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Find the max of scalar a the PEs in pelist
   result is also automatically broadcast to all PEs </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>a&nbsp;&nbsp;&nbsp;</tt></td><td> <tt>real</tt>   or <tt>integer</tt>, of 4-byte of 8-byte kind. </td>
</tr>
<tr>
<td valign="top" align="left"><tt>pelist&nbsp;&nbsp;&nbsp;</tt></td><td>   If <tt>pelist</tt>   is omitted, the context is assumed to be the
   current pelist. This call implies synchronization across the PEs in <tt>pelist</tt>, or the current pelist if <tt>pelist</tt>   is absent. </td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="mpp_sum"></a>
<h4>mpp_sum</h4>
<pre>
<b>call mpp_sum </b>( a, length, pelist )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd> 
<tt>MPP_TYPE_</tt>   corresponds to any 4-byte and 8-byte variant of <tt>integer, real, complex</tt>   variables, of rank 0 or 1. A
   contiguous block from a multi-dimensional array may be passed by its
   starting address and its length, as in <tt>f77</tt>.
   <br>
<br>
   Library reduction operators are not required or guaranteed to be
   bit-reproducible. In any case, changing the processor count changes
   the data layout, and thus very likely the order of operations. For
   bit-reproducible sums of distributed arrays, consider using the <tt>mpp_global_sum</tt>   routine provided by the <a href="models/bgrid_solo/fms_src/shared/mpp/mpp_domains.html"></a>   module.
   <br>
<br>
   The <tt>bit_reproducible</tt>   flag provided in earlier versions of
   this routine has been removed.
   <br>
<br>
   If <tt>pelist</tt>   is omitted, the context is assumed to be the
   current pelist. This call implies synchronization across the PEs in <tt>pelist</tt>, or the current pelist if <tt>pelist</tt>   is absent. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>length&nbsp;&nbsp;&nbsp;</tt></td><td></td>
</tr>
<tr>
<td valign="top" align="left"><tt>pelist&nbsp;&nbsp;&nbsp;</tt></td><td></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>INPUT/OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>a&nbsp;&nbsp;&nbsp;</tt></td><td></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="mpp_transmit"></a>
<h4>mpp_transmit</h4>
<pre>
<b>call mpp_transmit </b>( put_data, put_len, put_pe, get_data, get_len, get_pe )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd> 
<tt>MPP_TYPE_</tt>   corresponds to any 4-byte and 8-byte variant of <tt>integer, real, complex, logical</tt>   variables, of rank 0 or 1. A
   contiguous block from a multi-dimensional array may be passed by its
   starting address and its length, as in <tt>f77</tt>.
   <br>
<br> 
<tt>mpp_transmit</tt>   is currently implemented as asynchronous
   outward transmission and synchronous inward transmission. This follows
   the behaviour of <tt>shmem_put</tt>   and <tt>shmem_get</tt>. In MPI, it
   is implemented as <tt>mpi_isend</tt>   and <tt>mpi_recv</tt>. For most
   applications, transmissions occur in pairs, and are here accomplished
   in a single call.
   <br>
<br>
   The special PE designations <tt>NULL_PE</tt>, <tt>ANY_PE</tt>   and <tt>ALL_PES</tt>   are provided by use
   association.
   <br>
<br> 
<tt>NULL_PE</tt>: is used to disable one of the pair of
   transmissions.<br> 
<tt>ANY_PE</tt>: is used for unspecific remote
   destination. (Please note that <tt>put_pe=ANY_PE</tt>   has no meaning
   in the MPI context, though it is available in the SHMEM invocation. If
   portability is a concern, it is best avoided).<br> 
<tt>ALL_PES</tt>: is used for broadcast operations.
   <br>
<br>
   It is recommended that <a href="#mpp_broadcast"></a>   be used for
   broadcasts.
   <br>
<br>
   The following example illustrates the use of <tt>NULL_PE</tt>   and <tt>ALL_PES</tt>:
   <br>
<br> 
<pre>    real, dimension(n) :: a
    if( pe.EQ.0 )then
        do p = 1,npes-1
           call mpp_transmit( a, n, p, a, n, NULL_PE )
        end do
    else
        call mpp_transmit( a, n, NULL_PE, a, n, 0 )
    end if
    
    call mpp_transmit( a, n, ALL_PES, a, n, 0 )</pre>   The do loop and the broadcast operation above are equivalent.
   <br>
<br>
   Two overloaded calls <tt>mpp_send</tt>   and <tt>mpp_recv</tt>   have also been
   provided. <tt>mpp_send</tt>   calls <tt>mpp_transmit</tt>   with <tt>get_pe=NULL_PE</tt>. <tt>mpp_recv</tt>   calls <tt>mpp_transmit</tt>   with <tt>put_pe=NULL_PE</tt>. Thus
   the do loop above could be written more succinctly:
   <br>
<br> 
<pre>    if( pe.EQ.0 )then
        do p = 1,npes-1
           call mpp_send( a, n, p )
        end do
    else
        call mpp_recv( a, n, 0 )
    end if</pre> 
</dd>
<br>
<br>
</dl>
</li>
<li>
<a name="mpp_broadcast"></a>
<h4>mpp_broadcast</h4>
<pre>
<b>call mpp_broadcast </b>( data, length, from_pe, pelist )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   The <tt>mpp_broadcast</tt>   call has been added because the original
   syntax (using <tt>ALL_PES</tt>   in <tt>mpp_transmit</tt>) did not
   support a broadcast across a pelist.
   <br>
<br> 
<tt>MPP_TYPE_</tt>   corresponds to any 4-byte and 8-byte variant of <tt>integer, real, complex, logical</tt>   variables, of rank 0 or 1. A
   contiguous block from a multi-dimensional array may be passed by its
   starting address and its length, as in <tt>f77</tt>.
   <br>
<br>
   Global broadcasts through the <tt>ALL_PES</tt>   argument to <a href="#mpp_transmit"></a>   are still provided for
   backward-compatibility.
   <br>
<br>
   If <tt>pelist</tt>   is omitted, the context is assumed to be the
   current pelist. <tt>from_pe</tt>   must belong to the current
   pelist. This call implies synchronization across the PEs in <tt>pelist</tt>, or the current pelist if <tt>pelist</tt>   is absent. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>length&nbsp;&nbsp;&nbsp;</tt></td><td> </td>
</tr>
<tr>
<td valign="top" align="left"><tt>from_pe&nbsp;&nbsp;&nbsp;</tt></td><td> </td>
</tr>
<tr>
<td valign="top" align="left"><tt>pelist&nbsp;&nbsp;&nbsp;</tt></td><td> </td>
</tr>
</table>
</dd>
<br>
<dt>
<b>INPUT/OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>data(*)&nbsp;&nbsp;&nbsp;</tt></td><td> </td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="mpp_chksum"></a>
<h4>mpp_chksum</h4>
<pre> 
<b>mpp_chksum</b> ( var, pelist )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd> 
<tt>mpp_chksum</tt>   is a parallel checksum routine that returns an
   identical answer for the same array irrespective of how it has been
   partitioned across processors. <tt>LONG_KIND</tt>is the <tt>KIND</tt>   parameter corresponding to long integers (see discussion on
   OS-dependent preprocessor directives) defined in
   the header file <tt>os.h</tt>. <tt>MPP_TYPE_</tt>   corresponds to any
   4-byte and 8-byte variant of <tt>integer, real, complex, logical</tt>   variables, of rank 0 to 5.
   <br>
<br>
   Integer checksums on FP data use the F90 <tt>TRANSFER()</tt>   intrinsic.
   <br>
<br>
   The <a href="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/shared/chksum/chksum.html">serial checksum module</a>   is superseded
   by this function, and is no longer being actively maintained. This
   provides identical results on a single-processor job, and to perform
   serial checksums on a single processor of a parallel job, you only
   need to use the optional <tt>pelist</tt>   argument. <pre>     use mpp_mod
     integer :: pe, chksum
     real :: a(:)
     pe = mpp_pe()
     chksum = mpp_chksum( a, (/pe/) )</pre>   The additional functionality of <tt>mpp_chksum</tt>   over
   serial checksums is to compute the checksum across the PEs in <tt>pelist</tt>. The answer is guaranteed to be the same for
   the same distributed array irrespective of how it has been
   partitioned.
   <br>
<br>
   If <tt>pelist</tt>   is omitted, the context is assumed to be the
   current pelist. This call implies synchronization across the PEs in <tt>pelist</tt>, or the current pelist if <tt>pelist</tt>   is absent. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>pelist&nbsp;&nbsp;&nbsp;</tt></td><td> </td>
</tr>
<tr>
<td valign="top" align="left"><tt>var&nbsp;&nbsp;&nbsp;</tt></td><td> </td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="mpp_error"></a>
<h4>mpp_error</h4>
<pre>
<b>call mpp_error </b>( errortype, routine, errormsg )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   It is strongly recommended that all error exits pass through <tt>mpp_error</tt>   to assure the program fails cleanly. An individual
   PE encountering a <tt>STOP</tt>   statement, for instance, can cause the
   program to hang. The use of the <tt>STOP</tt>   statement is strongly
   discouraged.
   <br>
<br>
   Calling mpp_error with no arguments produces an immediate error
   exit, i.e: <pre>    call mpp_error
    call mpp_error(FATAL)</pre>   are equivalent.
   <br>
<br>
   The argument order <pre>    call mpp_error( routine, errormsg, errortype )</pre>   is also provided to support legacy code. In this version of the
   call, none of the arguments may be omitted.
   <br>
<br>
   The behaviour of <tt>mpp_error</tt>   for a <tt>WARNING</tt>   can be
   controlled with an additional call <tt>mpp_set_warn_level</tt>. <pre>    call mpp_set_warn_level(ERROR)</pre>   causes <tt>mpp_error</tt>   to treat <tt>WARNING</tt>   exactly like <tt>FATAL</tt>. <pre>    call mpp_set_warn_level(WARNING)</pre>   resets to the default behaviour described above.
   <br>
<br> 
<tt>mpp_error</tt>   also has an internal error state which
   maintains knowledge of whether a warning has been issued. This can be
   used at startup in a subroutine that checks if the model has been
   properly configured. You can generate a series of warnings using <tt>mpp_error</tt>, and then check at the end if any warnings has been
   issued using the function <tt>mpp_error_state()</tt>. If the value of
   this is <tt>WARNING</tt>, at least one warning has been issued, and
   the user can take appropriate action:
   <br>
<br> 
<pre>    if( ... )call mpp_error( WARNING, '...' )
    if( ... )call mpp_error( WARNING, '...' )
    if( ... )call mpp_error( WARNING, '...' )
    ...
    if( mpp_error_state().EQ.WARNING )call mpp_error( FATAL, '...' )</pre> 
</dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>errortype&nbsp;&nbsp;&nbsp;</tt></td><td>   One of <tt>NOTE</tt>, <tt>WARNING</tt>   or <tt>FATAL</tt>   (these definitions are acquired by use association). <tt>NOTE</tt>   writes <tt>errormsg</tt>   to <tt>STDOUT</tt>. <tt>WARNING</tt>   writes <tt>errormsg</tt>   to <tt>STDERR</tt>. <tt>FATAL</tt>   writes <tt>errormsg</tt>   to <tt>STDERR</tt>,
   and induces a clean error exit with a call stack traceback. </td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="mpp_init"></a>
<h4>mpp_init</h4>
<pre>
<b>call mpp_init </b>( flags )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Called to initialize the <tt>mpp_mod</tt>   package. It is recommended
   that this call be the first executed line in your program. It sets the
   number of PEs assigned to this run (acquired from the command line, or
   through the environment variable <tt>NPES</tt>), and associates an ID
   number to each PE. These can be accessed by calling <a href="#mpp_npes"></a>   and <a href="#mpp_pe"></a>. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>flags&nbsp;&nbsp;&nbsp;</tt></td><td> <tt>flags</tt>   can be set to <tt>MPP_VERBOSE</tt>   to
   have <tt>mpp_mod</tt>   keep you informed of what it's up to. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="stdin"></a>
<h4>stdin</h4>
<pre> 
<b>stdin</b> ()</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   This function, as well as stdout(), stderr(), stdlog(), returns the current 
   standard fortran unit numbers for
   input, output, error messages and log messages. Log messages, by
   convention, are written to the file <tt>logfile.out</tt>. </dd>
<br>
<br>
</dl>
</li>
<li>
<a name="mpp_exit"></a>
<h4>mpp_exit</h4>
<pre>
<b>call mpp_exit </b>()</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Called at the end of the run, or to re-initialize <tt>mpp_mod</tt>,
   should you require that for some odd reason.
   <br>
<br>
   This call implies synchronization across all PEs. </dd>
<br>
<br>
</dl>
</li>
<li>
<a name="mpp_pe"></a>
<h4>mpp_pe</h4>
<pre> 
<b>mpp_pe</b> ()</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   This returns the unique ID associated with a PE. This number runs
   between 0 and <tt>npes-1</tt>, where <tt>npes</tt>   is the total
   processor count, returned by <tt>mpp_npes</tt>. For a uniprocessor
   application this will always return 0. </dd>
<br>
<br>
</dl>
</li>
<li>
<a name="mpp_npes"></a>
<h4>mpp_npes</h4>
<pre> 
<b>mpp_npes</b> ()</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   This returns the number of PEs in the current pelist. For a
   uniprocessor application, this will always return 1. </dd>
<br>
<br>
</dl>
</li>
<li>
<a name="mpp_declare_pelist"></a>
<h4>mpp_declare_pelist</h4>
<pre>
<b>call mpp_declare_pelist </b>( pelist,name )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   This call is written specifically to accommodate a MPI restriction
   that requires a parent communicator to create a child communicator, In
   other words: a pelist cannot go off and declare a communicator, but
   every PE in the parent, including those not in pelist(:), must get
   together for the <tt>MPI_COMM_CREATE</tt>   call. The parent is
   typically <tt>MPI_COMM_WORLD</tt>, though it could also be a subset
   that includes all PEs in <tt>pelist</tt>.
   <br>
<br>
   The restriction does not apply to SMA but to have uniform code, you
   may as well call it.
   <br>
<br>
   This call implies synchronization across the PEs in the current
   pelist, of which <tt>pelist</tt>   is a subset. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>pelist&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, dimension(:)]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="mpp_set_current_pelist"></a>
<h4>mpp_set_current_pelist</h4>
<pre>
<b>call mpp_set_current_pelist </b>( pelist )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   This call sets the value of the current pelist, which is the
   context for all subsequent "global" calls where the optional <tt>pelist</tt>   argument is omitted. All the PEs that are to be in the
   current pelist must call it.
   <br>
<br>
   In MPI, this call may hang unless <tt>pelist</tt>   has been previous
   declared using <a href="#mpp_declare_pelist"></a>.
   <br>
<br>
   If the argument <tt>pelist</tt>   is absent, the current pelist is
   set to the "world" pelist, of all PEs in the job. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>pliest&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="mpp_clock_set_grain"></a>
<h4>mpp_clock_set_grain</h4>
<pre>
<b>call mpp_clock_set_grain </b>( grain )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   This routine and three other routines, mpp_clock_id, mpp_clock_begin(id),
   and mpp_clock_end(id) may be used to time parallel code sections, and
   extract parallel statistics. Clocks are identified by names, which
   should be unique in the first 32 characters. The <tt>mpp_clock_id</tt>   call initializes a clock of a given name and returns an integer <tt>id</tt>. This <tt>id</tt>   can be used by subsequent <tt>mpp_clock_begin</tt>   and <tt>mpp_clock_end</tt>   calls set around a
   code section to be timed. Example: <pre>    integer :: id
    id = mpp_clock_id( 'Atmosphere' )
    call mpp_clock_begin(id)
    call atmos_model()
    call mpp_clock_end()</pre>   Two flags may be used to alter the behaviour of <tt>mpp_clock</tt>. If the flag <tt>MPP_CLOCK_SYNC</tt>   is turned on
   by <tt>mpp_clock_id</tt>, the clock calls <tt>mpp_sync</tt>   across all
   the PEs in the current pelist at the top of the timed code section,
   but allows each PE to complete the code section (and reach <tt>mpp_clock_end</tt>) at different times. This allows us to measure
   load imbalance for a given code section. Statistics are written to <tt>stdout</tt>   by <tt>mpp_exit</tt>.
   <br>
<br>
   The flag <tt>MPP_CLOCK_DETAILED</tt>   may be turned on by <tt>mpp_clock_id</tt>   to get detailed communication
   profiles. Communication events of the types <tt>SEND, RECV, BROADCAST,
   REDUCE</tt>   and <tt>WAIT</tt>   are separately measured for data volume
   and time. Statistics are written to <tt>stdout</tt>   by <tt>mpp_exit</tt>, and individual PE info is also written to the file <tt>mpp_clock.out.####</tt>   where <tt>####</tt>   is the PE id given by <tt>mpp_pe</tt>.
   <br>
<br>
   The flags <tt>MPP_CLOCK_SYNC</tt>   and <tt>MPP_CLOCK_DETAILED</tt>   are
   integer parameters available by use association, and may be summed to
   turn them both on.
   <br>
<br>
   While the nesting of clocks is allowed, please note that turning on
   the non-optional flags on inner clocks has certain subtle issues.
   Turning on <tt>MPP_CLOCK_SYNC</tt>   on an inner
   clock may distort outer clock measurements of load imbalance. Turning
   on <tt>MPP_CLOCK_DETAILED</tt>   will stop detailed measurements on its
   outer clock, since only one detailed clock may be active at one time.
   Also, detailed clocks only time a certain number of events per clock
   (currently 40000) to conserve memory. If this array overflows, a
   warning message is printed, and subsequent events for this clock are
   not timed.
   <br>
<br>
   Timings are done using the <tt>f90</tt>   standard <tt>SYSTEM_CLOCK</tt>   intrinsic.
   <br>
<br>
   The resolution of SYSTEM_CLOCK is often too coarse for use except
   across large swaths of code. On SGI systems this is transparently
   overloaded with a higher resolution clock made available in a
   non-portable fortran interface made available by <tt>nsclock.c</tt>. This approach will eventually be extended to other
   platforms.
   <br>
<br>
   New behaviour added at the Havana release allows the user to embed
   profiling calls at varying levels of granularity all over the code,
   and for any particular run, set a threshold of granularity so that
   finer-grained clocks become dormant.
   <br>
<br>
   The threshold granularity is held in the private module variable <tt>clock_grain</tt>. This value may be modified by the call <tt>mpp_clock_set_grain</tt>, and affect clocks initiated by
   subsequent calls to <tt>mpp_clock_id</tt>. The value of <tt>clock_grain</tt>   is set to an arbitrarily large number initially.
   <br>
<br>
   Clocks initialized by <tt>mpp_clock_id</tt>   can set a new optional
   argument <tt>grain</tt>   setting their granularity level. Clocks check
   this level against the current value of <tt>clock_grain</tt>, and are
   only triggered if they are <i>at or below ("coarser than")</i>   the
   threshold. Finer-grained clocks are dormant for that run.
   <br>
<br>
   Note that subsequent changes to <tt>clock_grain</tt>   do not
   change the status of already initiated clocks, and that if the
   optional <tt>grain</tt>   argument is absent, the clock is always
   triggered. This guarantees backward compatibility. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>grain&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="mpp_sync"></a>
<h4>mpp_sync</h4>
<pre>
<b>call mpp_sync </b>( pelist )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Synchronizes PEs at this point in the execution. If <tt>pelist</tt>   is omitted all PEs are synchronized. This can
   be expensive on many systems, and should be avoided if possible. Under
   MPI, we do not call <tt>MPI_BARRIER</tt>, as you might
   expect. This is because this call can be prohibitively slow on many
   systems. Instead, we perform the same operation as <tt>mpp_sync_self</tt>, i.e all participating PEs wait for
   completion of all their outstanding non-blocking operations.
   <br>
<br>
   If <tt>pelist</tt>   is omitted, the context is assumed to be the
   current pelist. This call implies synchronization across the PEs in <tt>pelist</tt>, or the current pelist if <tt>pelist</tt>   is absent. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>pelist&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, dimension(:)]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="mpp_sync_self"></a>
<h4>mpp_sync_self</h4>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd> 
<tt>mpp_transmit</tt>   is implemented as asynchronous <tt>put/send</tt>   and synchronous <tt>get/recv</tt>. <tt>mpp_sync_self</tt>   guarantees that outstanding
   asynchronous operations from the calling PE are complete. If <tt>pelist</tt>   is supplied, <tt>mpp_sync_self</tt>   checks only for
   outstanding puts to the PEs in <tt>pelist</tt>.
   <br>
<br>
   If <tt>pelist</tt>   is omitted, the context is assumed to be the
   current pelist. This call implies synchronization across the PEs in <tt>pelist</tt>, or the current pelist if <tt>pelist</tt>   is absent. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>pelist&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, dimension(:)]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="mpp_malloc"></a>
<h4>mpp_malloc</h4>
<pre>
<b>call mpp_malloc </b>( ptr, newlen, len )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   This routine is used on SGI systems when <tt>mpp_mod</tt>   is
   invoked in the SHMEM library. It ensures that dynamically allocated
   memory can be used with <tt>shmem_get</tt>   and <tt>shmem_put</tt>. This is called <i>symmetric
   allocation</i>   and is described in the <tt>intro_shmem</tt>   man page. <tt>ptr</tt>   is a <i>Cray
   pointer</i>   (see the section on <a href="#PORTABILITY">portability</a>).  The operation can be expensive
   (since it requires a global barrier). We therefore attempt to re-use
   existing allocation whenever possible. Therefore <tt>len</tt>   and <tt>ptr</tt>   must have the <tt>SAVE</tt>   attribute
   in the calling routine, and retain the information about the last call
   to <tt>mpp_malloc</tt>. Additional memory is symmetrically
   allocated if and only if <tt>newlen</tt>   exceeds <tt>len</tt>.
   <br>
<br>
   This is never required on Cray PVP or MPP systems. While the T3E
   manpages do talk about symmetric allocation, <tt>mpp_mod</tt>   is coded to remove this restriction.
   <br>
<br>
   It is never required if <tt>mpp_mod</tt>   is invoked in MPI.
   <br>
<br>
   This call implies synchronization across all PEs. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>ptr&nbsp;&nbsp;&nbsp;</tt></td><td>   a cray pointer, points to a dummy argument in this routine. </td>
</tr>
<tr>
<td valign="top" align="left"><tt>newlen&nbsp;&nbsp;&nbsp;</tt></td><td>   the required allocation length for the pointer ptr <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>len&nbsp;&nbsp;&nbsp;</tt></td><td>   the current allocation (0 if unallocated). <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="mpp_set_stack_size"></a>
<h4>mpp_set_stack_size</h4>
<pre>
<b>call mpp_set_stack_size </b>(n)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd> 
<tt>mpp_mod</tt>   maintains a private internal array called <tt>mpp_stack</tt>   for private workspace. This call sets the length,
   in words, of this array. 
   <br>
<br>
   The <tt>mpp_init</tt>   call sets this
   workspace length to a default of 32768, and this call may be used if a
   longer workspace is needed.
   <br>
<br>
   This call implies synchronization across all PEs.
   <br>
<br>
   This workspace is symmetrically allocated, as required for
   efficient communication on SGI and Cray MPP systems. Since symmetric
   allocation must be performed by <i>all</i>   PEs in a job, this call
   must also be called by all PEs, using the same value of <tt>n</tt>. Calling <tt>mpp_set_stack_size</tt>   from a subset of PEs,
   or with unequal argument <tt>n</tt>, may cause the program to hang.
   <br>
<br>
   If any MPP call using <tt>mpp_stack</tt>   overflows the declared
   stack array, the program will abort with a message specifying the
   stack length that is required. Many users wonder why, if the required
   stack length can be computed, it cannot also be specified at that
   point. This cannot be automated because there is no way for the
   program to know if all PEs are present at that call, and with equal
   values of <tt>n</tt>. The program must be rerun by the user with the
   correct argument to <tt>mpp_set_stack_size</tt>, called at an
   appropriate point in the code where all PEs are known to be present. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>n&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
</table>
</dd>
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
<div>None.<br>
<br>
</div>
<!-- END ERROR MESSAGES -->
<a name="REFERENCES"></a>
<hr>
<h4>REFERENCES</h4>
<!-- BEGIN REFERENCES -->
<div>
        None.
      </div>
<br>
<!-- END REFERENCES -->
<a name="COMPILER SPECIFICS"></a>
<hr>
<h4>COMPILER SPECIFICS</h4>
<!-- BEGIN COMPILER SPECIFICS -->
<div>
<dl>
<dt></dt>
<dd>   Any module or program unit using <tt>mpp_mod</tt>   must contain the line <pre>    use mpp_mod</pre>   The source file for <tt>mpp_mod</tt>   is <a href="ftp://ftp.gfdl.gov/pub/vb/mpp/mpp.F90"></a>.
   Activate the preprocessor flag <tt>-Duse_libSMA</tt>   to invoke
   the SHMEM library, or <tt>-Duse_libMPI</tt>   to invoke the MPI
   library. Global translation of preprocessor macros is required. This
   required the activation of the <tt>-F</tt>   flag on Cray systems
   and the <tt>-ftpp -macro_expand</tt>   flags on SGI systems. On
   non-SGI/Cray systems, please consult the f90 manpage for the
   equivalent flag.
   <br>
<br>
   On Cray PVP systems, <i>all</i>   routines in a message-passing
   program must be compiled with <tt>-a taskcommon</tt>.
   <br>
<br>
   On SGI systems, it is required to use 4-byte integers and 8-byte
   reals, and the 64-bit ABI (<tt>-i4 -r8 -64 -mips4</tt>). It is also
   required on SGI systems to link the following libraries explicitly:
   one of <tt>-lmpi</tt>   and <tt>-lsma</tt>, depending on
   whether you wish to use the SHMEM or MPI implementations; and <tt>-lexc</tt>). On Cray systems, all the required flags are
   default.
   <br>
<br>
   On SGI, use MIPSPro f90 7.3.1.2 or higher.
   <br>
<br>
   On Cray, use cf90 3.0.0.0 or higher.
   <br>
<br>
   On either, use the message-passing toolkit MPT 1.2 or higher.
   <br>
<br>
   The declaration <tt>MPI_INTEGER8</tt>   for 8-byte integers
   was provided by <tt>mpp_mod</tt>   because it was absent in early
   releases of the Message Passing Toolkit. It has since been included
   there, and the declaration in <tt>mpp_mod</tt>   commented
   out. This declaration may need to be reinstated if you get a compiler
   error from this (i.e you are using a superseded version of the MPT).
   <br>
<br>
   By turning on the cpp flag <tt>-Dtest_mpp</tt>   and compiling <tt>mpp_mod</tt>   by itself, you may create a test program to
   exercise certain aspects of <tt>mpp_mod</tt>, e.g
   <br>
<br> 
<pre>    f90 -F -Duse_libSMA -Dtest_mpp mpp.F90
    mpprun -n4 a.out</pre>   runs a 4-PE test on a t3e. </dd>
</dl>
</div>
<br>
<!-- END COMPILER SPECIFICS -->
<a name="PRECOMPILER OPTIONS"></a>
<hr>
<h4>PRECOMPILER OPTIONS</h4>
<!-- BEGIN PRECOMPILER OPTIONS -->
<div>
<dl>
<dt></dt>
<dd>   While the SHMEM library is currently available only on SGI/Cray
   systems, <tt>mpp_mod</tt>   can be used on any other system with
   a standard-compliant f90 compiler and MPI library. SHMEM is now
   becoming available on other systems as well.
   <br>
<br>
   There are some <a href="">OS-dependent
   pre-processor directives</a>   that you might need to modify on
   non-SGI/Cray systems and compilers.
   <br>
<br>
   On SGI systems, the <tt>f90</tt>   standard <tt>SYSTEM_CLOCK</tt>   intrinsic is overloaded with a non-portable fortran interface to a
   higher-precision clock. This is distributed with the MPP package as <tt>nsclock.c</tt>. This approach will eventually be extended to other
   platforms, since the resolution of the default clock is often too
   coarse for our needs. </dd>
</dl>
</div>
<br>
<!-- END PRECOMPILER OPTIONS -->
<a name="LOADER OPTIONS"></a>
<hr>
<h4>LOADER OPTIONS</h4>
<!-- BEGIN LOADER -->
<div>
<p>   The    source consists of the main source file    and also requires the following include files:
       (when compiled with )    (when compiled with )         GFDL users can check it out of the main CVS repository as part of
   the    CVS module. The current public tag is .
   External users can download the latest    package . Public access
   to the GFDL CVS repository will soon be made available. </p>
<pre>        
</pre>
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
<p>   The <tt>SYSTEM_CLOCK</tt>   intrinsic has a limited range before the
   clock rolls over. The maximum time interval that may be measured
   before rollover depends on the default integer precision, and is <tt>COUNT_MAX/COUNT_RATE</tt>   seconds. Timing a code section longer
   than this interval will give incorrect results. The <tt>mpp</tt>   entry in the logfile reports the rollover time interval. Note that
   this is a limitation, or "feature" of the <tt>f90 SYSTEM_CLOCK</tt>   intrinsic. </p>
</div>
<br>
<!-- END KNOWN BUGS -->
<a name="NOTES"></a>
<hr>
<h4>NOTES</h4>
<!-- BEGIN NOTES -->
<div>None.<br>
</div>
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
