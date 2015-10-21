//*****************************************************************************
//
// Project:	PROJ.4
// Purpose:	Functions for handling individual PJ_GRIDINFO's. Includes
//			loaders for all formats but CTABLE (in nad_init.cs).
// Author:	Frank Warmerdam, warmerdam@pobox.com
//
//*****************************************************************************
// Copyright (c) 2000, Frank Warmerdam <warmerdam@pobox.com>
// Copyright (c) 2009-2015 by the Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//*****************************************************************************

using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Gridshift
{
	public class PJ_GRIDINFO
	{
		public string gridname;		// identifying name of grid, eg "conus" or ntv2_0.gsb
		public string filename;		// full path to filename

		public string format;		// format of this grid, ie "ctable", "ntv1", "ntv2" or "missing".

		public long grid_offset;	// offset in file, for delayed loading

		public CTABLE ct;

		public PJ_GRIDINFO next;
		public PJ_GRIDINFO child;
	}

	public static partial class Grid
	{
		static object gridlock=new object();

		//**********************************************************************
		//							swap_words()
		//
		//		Convert the byte order of the given word(s) in place.
		//**********************************************************************
		//static int byte_order_test=1;
		static bool IS_LSB=true;//	(1 == ((unsigned char *) (&byte_order_test))[0])

		static void swap_words(byte[] data, int word_size, int word_count)
		{
			swap_words(data, 0, word_size, word_count);
		}

		static void swap_words(byte[] data, int offset, int word_size, int word_count)
		{
			for(int word=0; word<word_count; word++)
			{
				for(int i=0; i<word_size/2; i++)
				{
					byte t=data[offset+i];
					data[i]=data[offset+word_size-i-1];
					data[offset+word_size-i-1]=t;
				}

				offset+=word_size;
			}
		}

		//**********************************************************************
		//							pj_gridinfo_free()
		//**********************************************************************
		public static void pj_gridinfo_free(projCtx ctx, PJ_GRIDINFO gi)
		{
			if(gi==null) return;

			if(gi.child!=null)
			{
				PJ_GRIDINFO child, next;

				for(child=gi.child; child!=null; child=next)
				{
					next=child.next;
					pj_gridinfo_free(ctx, child);
					child=null;
				}
			}

			if(gi.ct!=null) nad_free(gi.ct);
			gi.ct=null;
		}

		//**********************************************************************
		//							pj_gridinfo_load()
		//
		//		This function is intended to implement delayed loading of
		//		the data contents of a grid file. The header and related
		//		stuff are loaded by pj_gridinfo_init().
		//**********************************************************************
		public static bool pj_gridinfo_load(projCtx ctx, PJ_GRIDINFO gi)
		{
			if(gi==null||gi.ct==null) return false;

			lock(gridlock)
			{
				if(gi.ct.cvs!=null) return true;

				CTABLE ct_tmp=gi.ct.Clone();

				// --------------------------------------------------------------------
				//		Original platform specific CTable format.
				// --------------------------------------------------------------------
				if(gi.format=="ctable")
				{
					FileStream fid=Proj.pj_open_lib(ctx, gi.filename, FileAccess.Read);
					if(fid==null)
					{
						Proj.pj_ctx_set_errno(ctx, -38);
						return false;
					}

					bool result=nad_ctable_load(ctx, ct_tmp, fid);

					fid.Close();

					gi.ct.cvs=ct_tmp.cvs;

					return result;
				}
				// --------------------------------------------------------------------
				//		CTable2 format.
				// --------------------------------------------------------------------
				else if(gi.format=="ctable2")
				{
					FileStream fid=Proj.pj_open_lib(ctx, gi.filename, FileAccess.Read);
					if(fid==null)
					{
						Proj.pj_ctx_set_errno(ctx, -38);
						return false;
					}

					bool result=nad_ctable2_load(ctx, ct_tmp, fid);

					fid.Close();

					gi.ct.cvs=ct_tmp.cvs;

					return result;
				}

				// --------------------------------------------------------------------
				//		NTv1 format.
				//		We process one line at a time. Note that the array storage
				//		direction (e-w) is different in the NTv1 file and what
				//		the CTABLE is supposed to have. The phi/lam are also
				//		reversed, and we have to be aware of byte swapping.
				// --------------------------------------------------------------------
				if(gi.format=="ntv1")
				{
					FileStream fid=Proj.pj_open_lib(ctx, gi.filename, FileAccess.Read);

					if(fid==null)
					{
						Proj.pj_ctx_set_errno(ctx, -38);
						return false;
					}

					fid.Seek(gi.grid_offset, SeekOrigin.Begin);

					byte[] row_buf;
					try
					{
						row_buf=new byte[gi.ct.lim.lam*2*sizeof(double)];
						ct_tmp.cvs=new LP[gi.ct.lim.lam*gi.ct.lim.phi];
					}
					catch
					{
						row_buf=null;
						ct_tmp.cvs=null;

						Proj.pj_ctx_set_errno(ctx, -38);
						return false;
					}

					for(int row=0; row<gi.ct.lim.phi; row++)
					{
						try
						{
							if(fid.Read(row_buf, 0, gi.ct.lim.lam*2*sizeof(double))!=gi.ct.lim.lam*2*sizeof(double))
							{
								row_buf=null;
								ct_tmp.cvs=null;

								Proj.pj_ctx_set_errno(ctx, -38);
								return false;
							}
						}
						catch
						{
							row_buf=null;
							ct_tmp.cvs=null;

							Proj.pj_ctx_set_errno(ctx, -38);
							return false;
						}

						if(IS_LSB) swap_words(row_buf, 8, gi.ct.lim.lam*2);

						// convert seconds to radians
						int diff_seconds=0;

						for(int i=0; i<gi.ct.lim.lam; i++)
						{
							int cvs=row*gi.ct.lim.lam+(gi.ct.lim.lam-i-1);

							ct_tmp.cvs[cvs].phi=BitConverter.ToDouble(row_buf, (diff_seconds++)*sizeof(double))*((Proj.PI/180.0)/3600.0);
							ct_tmp.cvs[cvs].lam=BitConverter.ToDouble(row_buf, (diff_seconds++)*sizeof(double))*((Proj.PI/180.0)/3600.0);
						}
					}

					row_buf=null;

					fid.Close();

					gi.ct.cvs=ct_tmp.cvs;

					return true;
				}

				// --------------------------------------------------------------------
				//		NTv2 format.
				//		We process one line at a time. Note that the array storage
				//		direction (e-w) is different in the NTv2 file and what
				//		the CTABLE is supposed to have. The phi/lam are also
				//		reversed, and we have to be aware of byte swapping.
				// --------------------------------------------------------------------
				if(gi.format=="ntv2")
				{
#if DEBUG
					Console.Error.WriteLine("NTv2 - loading grid "+gi.ct.id);
#endif

					FileStream fid=Proj.pj_open_lib(ctx, gi.filename, FileAccess.Read);

					if(fid==null)
					{
						Proj.pj_ctx_set_errno(ctx, -38);
						return false;
					}

					fid.Seek(gi.grid_offset, SeekOrigin.Begin);

					byte[] row_buf;
					try
					{
						row_buf=new byte[gi.ct.lim.lam*4*sizeof(float)];
						ct_tmp.cvs=new LP[gi.ct.lim.lam*gi.ct.lim.phi];
					}
					catch
					{
						row_buf=null;
						ct_tmp.cvs=null;

						Proj.pj_ctx_set_errno(ctx, -38);
						return false;
					}

					for(int row=0; row<gi.ct.lim.phi; row++)
					{
						try
						{
							if(fid.Read(row_buf, 0, gi.ct.lim.lam*4*sizeof(float))!=gi.ct.lim.lam*4*sizeof(float))
							{
								row_buf=null;
								ct_tmp.cvs=null;

								Proj.pj_ctx_set_errno(ctx, -38);
								return false;
							}
						}
						catch
						{
							row_buf=null;
							ct_tmp.cvs=null;

							Proj.pj_ctx_set_errno(ctx, -38);
							return false;
						}

						if(!IS_LSB) swap_words(row_buf, 4, gi.ct.lim.lam*4);

						// convert seconds to radians
						int diff_seconds=0;

						for(int i=0; i<gi.ct.lim.lam; i++)
						{
							int cvs=row*gi.ct.lim.lam+(gi.ct.lim.lam-i-1);

							ct_tmp.cvs[cvs].phi=BitConverter.ToSingle(row_buf, (diff_seconds++)*sizeof(float))*((Proj.PI/180.0)/3600.0);
							ct_tmp.cvs[cvs].lam=BitConverter.ToSingle(row_buf, (diff_seconds++)*sizeof(float))*((Proj.PI/180.0)/3600.0);
							diff_seconds+=2; // skip accuracy values
						}
					}

					row_buf=null;

					fid.Close();

					gi.ct.cvs=ct_tmp.cvs;

					return true;
				}

				// --------------------------------------------------------------------
				//		GTX format.
				// --------------------------------------------------------------------
				if(gi.format=="gtx")
				{
					int words=gi.ct.lim.lam*gi.ct.lim.phi;
					FileStream fid=Proj.pj_open_lib(ctx, gi.filename, FileAccess.Read);
					if(fid==null)
					{
						Proj.pj_ctx_set_errno(ctx, -38);
						return false;
					}

					fid.Seek(gi.grid_offset, SeekOrigin.Begin);

					byte[] buf;

					try
					{
						buf=new byte[words*sizeof(float)];
						ct_tmp.cvs=new LP[words/2];
					}
					catch
					{
						Proj.pj_ctx_set_errno(ctx, -38);
						return false;
					}

					try
					{
						if(fid.Read(buf, 0, words*sizeof(float))!=words*sizeof(float))
						{
							buf=null;
							ct_tmp.cvs=null;
							return false;
						}
					}
					catch
					{
						buf=null;
						ct_tmp.cvs=null;
						return false;
					}

					if(IS_LSB) swap_words(buf, 4, words);

					for(int i=0; i<words; i+=2)
					{
						ct_tmp.cvs[i/2].phi=BitConverter.ToSingle(buf, i*sizeof(float));
						ct_tmp.cvs[i/2].lam=BitConverter.ToSingle(buf, (i+1)*sizeof(float));
					}

					fid.Close();

					gi.ct.cvs=ct_tmp.cvs;

					return true;
				}

				return false;
			} // lock(gridlock)
		}

		//**********************************************************************
		//						pj_gridinfo_parent()
		//
		//		Seek a parent grid file by name from a grid list.
		//**********************************************************************
		static PJ_GRIDINFO pj_gridinfo_parent(PJ_GRIDINFO gilist, string name)
		{
			while(gilist!=null)
			{
				if(gilist.ct.id==name) return gilist;

				if(gilist.child!=null)
				{
					PJ_GRIDINFO parent=pj_gridinfo_parent(gilist.child, name);
					if(parent!=null) return parent;
				}

				gilist=gilist.next;
			}

			return gilist;
		}

		//**********************************************************************
		//						pj_gridinfo_init_ntv2()
		//
		//		Load a ntv2 (.gsb) file.
		//**********************************************************************
		static bool pj_gridinfo_init_ntv2(projCtx ctx, FileStream fid, PJ_GRIDINFO gilist)
		{
			byte[] header=new byte[11*16];

			// --------------------------------------------------------------------
			//		Read the overview header.
			// --------------------------------------------------------------------
			try
			{
				if(fid.Read(header, 0, header.Length)!=header.Length)
				{
					header=null;
					Proj.pj_ctx_set_errno(ctx, -38);
					return false;
				}
			}
			catch
			{
				header=null;
				Proj.pj_ctx_set_errno(ctx, -38);
				return false;
			}

			// --------------------------------------------------------------------
			//		Byte swap interesting fields if needed.
			// --------------------------------------------------------------------
			if(!IS_LSB)
			{
				swap_words(header, 8, 4, 1);
				swap_words(header, 8+16, 4, 1);
				swap_words(header, 8+32, 4, 1);
				swap_words(header, 8+7*16, 8, 1);
				swap_words(header, 8+8*16, 8, 1);
				swap_words(header, 8+9*16, 8, 1);
				swap_words(header, 8+10*16, 8, 1);
			}

			// --------------------------------------------------------------------
			//		Get the subfile count out ... all we really use for now.
			// --------------------------------------------------------------------
			int num_subfiles=BitConverter.ToInt32(header, 8+32);

			// ====================================================================
			//		Step through the subfiles, creating a PJ_GRIDINFO for each.
			// ====================================================================
			for(int subfile=0; subfile<num_subfiles; subfile++)
			{
				// --------------------------------------------------------------------
				//		Read header.
				// --------------------------------------------------------------------
				try
				{
					if(fid.Read(header, 0, header.Length)!=header.Length)
					{
						header=null;
						Proj.pj_ctx_set_errno(ctx, -38);
						return false;
					}
				}
				catch
				{
					header=null;
					Proj.pj_ctx_set_errno(ctx, -38);
					return false;
				}

				if(Encoding.ASCII.GetString(header, 0, 8).StartsWith("SUB_NAME"))
				{
					Proj.pj_ctx_set_errno(ctx, -38);
					return false;
				}

				// --------------------------------------------------------------------
				//		Byte swap interesting fields if needed.
				// --------------------------------------------------------------------
				if(!IS_LSB)
				{
					swap_words(header, 8+16*4, 8, 1);
					swap_words(header, 8+16*5, 8, 1);
					swap_words(header, 8+16*6, 8, 1);
					swap_words(header, 8+16*7, 8, 1);
					swap_words(header, +8+16*8, 8, 1);
					swap_words(header, 8+16*9, 8, 1);
					swap_words(header, 8+16*10, 4, 1);
				}

				// --------------------------------------------------------------------
				//		Initialize a corresponding "ct" structure.
				// --------------------------------------------------------------------
				CTABLE ct=new CTABLE();
				ct.id=Encoding.ASCII.GetString(header, 8, 8);

				ct.ll.lam=-BitConverter.ToDouble(header, 7*16+8);	// W_LONG
				ct.ll.phi=BitConverter.ToDouble(header, 4*16+8);	// S_LAT

				LP ur;
				ur.lam=-BitConverter.ToDouble(header, 6*16+8);		// E_LONG
				ur.phi=BitConverter.ToDouble(header, 5*16+8);		// N_LAT

				ct.del.lam=BitConverter.ToDouble(header, 9*16+8);
				ct.del.phi=BitConverter.ToDouble(header, 8*16+8);

				ct.lim.lam=(int)(Math.Abs(ur.lam-ct.ll.lam)/ct.del.lam+0.5)+1;
				ct.lim.phi=(int)(Math.Abs(ur.phi-ct.ll.phi)/ct.del.phi+0.5)+1;

#if DEBUG
				Console.Error.WriteLine("NTv2 {0} {1}x{2}: LL=({3},{4}) UR=({5},{6})", ct.id,
					ct.lim.lam, ct.lim.phi, ct.ll.lam/3600.0, ct.ll.phi/3600.0, ur.lam/3600.0, ur.phi/3600.0);
#endif

				ct.ll.lam*=Proj.DEG_TO_RAD/3600.0;
				ct.ll.phi*=Proj.DEG_TO_RAD/3600.0;
				ct.del.lam*=Proj.DEG_TO_RAD/3600.0;
				ct.del.phi*=Proj.DEG_TO_RAD/3600.0;

				int gs_count=BitConverter.ToInt32(header, 8+16*10);
				if(gs_count!=ct.lim.lam*ct.lim.phi)
				{
					Console.Error.WriteLine("GS_COUNT({0}) does not match expected cells ({1}x{2}={3})",
						gs_count, ct.lim.lam, ct.lim.phi, ct.lim.lam*ct.lim.phi);
					Proj.pj_ctx_set_errno(ctx, -38);
					return false;
				}

				ct.cvs=null;

				PJ_GRIDINFO gi;

				// --------------------------------------------------------------------
				//		Create a new gridinfo for this if we aren't processing the
				//		1st subfile, and initialize our grid info.
				// --------------------------------------------------------------------
				if(subfile==0) gi=gilist;
				else
				{
					gi=new PJ_GRIDINFO();
					gi.gridname=gilist.gridname;
					gi.filename=gilist.filename;
					gi.next=null;
				}

				gi.ct=ct;
				gi.format="ntv2";
				gi.grid_offset=fid.Position;

				// --------------------------------------------------------------------
				//		Attach to the correct list or sublist.
				// --------------------------------------------------------------------
				if(Encoding.ASCII.GetString(header, 24, 4)=="NONE")
				{
					if(gi!=gilist)
					{
						PJ_GRIDINFO lnk=gilist;
						while(lnk.next!=null) lnk=lnk.next;
						lnk.next=gi;
					}
				}
				else
				{
					PJ_GRIDINFO gp=pj_gridinfo_parent(gilist, Encoding.ASCII.GetString(header, 24, 8));

					if(gp==null)
					{
#if DEBUG
						Console.Error.WriteLine("pj_gridinfo_init_ntv2(): failed to find parent {0} for {1}.",
							Encoding.ASCII.GetString(header, 24, 8), gi.ct.id);
#endif

						PJ_GRIDINFO lnk=gilist;
						while(lnk.next!=null) lnk=lnk.next;
						lnk.next=gi;
					}
					else if(gp.child==null) gp.child=gi;
					else
					{
						PJ_GRIDINFO lnk=gp.child;
						while(lnk.next!=null) lnk=lnk.next;
						lnk.next=gi;
					}
				}

				// --------------------------------------------------------------------
				//		Seek past the data.
				// --------------------------------------------------------------------
				fid.Seek(gs_count*16, SeekOrigin.Current);
			}

			return true;
		}

		//**********************************************************************
		//						pj_gridinfo_init_ntv1()
		//
		//		Load an NTv1 style Canadian grid shift file.
		//**********************************************************************
		static bool pj_gridinfo_init_ntv1(projCtx ctx, FileStream fid, PJ_GRIDINFO gi)
		{
			byte[] header=new byte[176];

			// --------------------------------------------------------------------
			//		Read the header.
			// --------------------------------------------------------------------
			try
			{
				if(fid.Read(header, 0, header.Length)!=header.Length)
				{
					header=null;
					Proj.pj_ctx_set_errno(ctx, -38);
					return false;
				}
			}
			catch
			{
				header=null;
				Proj.pj_ctx_set_errno(ctx, -38);
				return false;
			}

			// --------------------------------------------------------------------
			//		Regularize fields of interest.
			// --------------------------------------------------------------------
			if(IS_LSB)
			{
				swap_words(header, 8, 4, 1);
				swap_words(header, 24, 8, 1);
				swap_words(header, 40, 8, 1);
				swap_words(header, 56, 8, 1);
				swap_words(header, 72, 8, 1);
				swap_words(header, 88, 8, 1);
				swap_words(header, 104, 8, 1);
			}

			if(BitConverter.ToInt32(header, 8)!=12)
			{
				Console.Error.WriteLine("NTv1 grid shift file has wrong record count, corrupt?");
				Proj.pj_ctx_set_errno(ctx, -38);
				return false;
			}

			// --------------------------------------------------------------------
			//		Fill in CTABLE structure.
			// --------------------------------------------------------------------
			CTABLE ct=new CTABLE();
			ct.id="NTv1 Grid Shift File";

			ct.ll.lam=-BitConverter.ToDouble(header, 72);
			ct.ll.phi=BitConverter.ToDouble(header, 24);

			LP ur;
			ur.lam=-BitConverter.ToDouble(header, 56);
			ur.phi=BitConverter.ToDouble(header, 40);

			ct.del.lam=BitConverter.ToDouble(header, 104);
			ct.del.phi=BitConverter.ToDouble(header, 88);

			ct.lim.lam=(int)(Math.Abs(ur.lam-ct.ll.lam)/ct.del.lam+0.5)+1;
			ct.lim.phi=(int)(Math.Abs(ur.phi-ct.ll.phi)/ct.del.phi+0.5)+1;

#if DEBUG
			Console.Error.WriteLine("NTv1 {0}x{1}: LL=({2},{3}) UR=({4},{5})", 
				ct.lim.lam, ct.lim.phi, ct.ll.lam, ct.ll.phi, ur.lam, ur.phi);
#endif

			ct.ll.lam*=Proj.DEG_TO_RAD;
			ct.ll.phi*=Proj.DEG_TO_RAD;
			ct.del.lam*=Proj.DEG_TO_RAD;
			ct.del.phi*=Proj.DEG_TO_RAD;
			ct.cvs=null;

			gi.ct=ct;
			gi.grid_offset=fid.Position;
			gi.format="ntv1";

			return true;
		}

		//**********************************************************************
		//						pj_gridinfo_init_gtx()
		//
		//		Load a NOAA .gtx vertical datum shift file.
		//**********************************************************************
		static bool pj_gridinfo_init_gtx(projCtx ctx, FileStream fid, PJ_GRIDINFO gi)
		{
			byte[] header=new byte[40];
			CTABLE ct;
			double xorigin, yorigin, xstep, ystep;
			int rows, columns;

			// --------------------------------------------------------------------
			//		Read the header.
			// --------------------------------------------------------------------
			if(fid.Read(header, 0, 40)!=40)
			{
				Proj.pj_ctx_set_errno(ctx, -38);
				return false;
			}

			// --------------------------------------------------------------------
			//		Regularize fields of interest and extract.
			// --------------------------------------------------------------------
			if(IS_LSB)
			{
				swap_words(header, 0, 8, 4);
				swap_words(header, 32, 4, 2);
			}

			yorigin=BitConverter.ToDouble(header, 0);
			xorigin=BitConverter.ToDouble(header, 8);
			ystep=BitConverter.ToDouble(header, 16);
			xstep=BitConverter.ToDouble(header, 24);

			rows=BitConverter.ToInt32(header, 32);
			columns=BitConverter.ToInt32(header, 36);

			if(xorigin<-360||xorigin>360
				||yorigin<-90||yorigin>90)
			{
				Console.WriteLine("gtx file header has invalid extents, corrupt?");
				Proj.pj_ctx_set_errno(ctx, -38);
				return false;
			}

			// --------------------------------------------------------------------
			//		Fill in CTABLE structure.
			// --------------------------------------------------------------------
			ct=new CTABLE();
			ct.id="GTX Vertical Grid Shift File";

			ct.ll.lam=xorigin;
			ct.ll.phi=yorigin;
			ct.del.lam=xstep;
			ct.del.phi=ystep;
			ct.lim.lam=columns;
			ct.lim.phi=rows;

			// some GTX files come in 0-360 and we shift them back into the
			// expected -180 to 180 range if possible. This does not solve
			// problems with grids spanning the dateline.
			if(ct.ll.lam>=180.0) ct.ll.lam-=360.0;

#if !DEBUG
			if(ct.ll.lam>=0.0&&ct.ll.lam+ct.del.lam*ct.lim.lam>180.0)
				Console.Error.WriteLine("This GTX spans the dateline! This will cause problems.");

			Console.Error.WriteLine("GTX {0}x{1}: LL=({2},{3}) UR=({4},{5})", ct.lim.lam, ct.lim.phi, ct.ll.lam, ct.ll.phi, ct.ll.lam+(columns-1)*xstep, ct.ll.phi+(rows-1)*ystep);
#endif

			ct.ll.lam*=Proj.DEG_TO_RAD;
			ct.ll.phi*=Proj.DEG_TO_RAD;
			ct.del.lam*=Proj.DEG_TO_RAD;
			ct.del.phi*=Proj.DEG_TO_RAD;
			ct.cvs=null;

			gi.ct=ct;
			gi.grid_offset=40;
			gi.format="gtx";

			return true;
		}

		//**********************************************************************
		//							pj_gridinfo_init()
		//
		//		Open and parse header details from a datum gridshift file
		//		returning a list of PJ_GRIDINFOs for the grids in that
		//		file. This superceeds use of nad_init() for modern
		//		applications.
		//**********************************************************************
		public static PJ_GRIDINFO pj_gridinfo_init(projCtx ctx, string gridname)
		{
			Libc.errno=Proj.pj_errno=0;
			ctx.last_errno=0;

			// --------------------------------------------------------------------
			//		Initialize a GRIDINFO with stub info we would use if it
			//		cannot be loaded.
			// --------------------------------------------------------------------
			PJ_GRIDINFO gilist=new PJ_GRIDINFO();
			gilist.gridname=gridname;
			gilist.filename=null;
			gilist.format="missing";
			gilist.grid_offset=0;
			gilist.ct=null;
			gilist.next=null;

			// --------------------------------------------------------------------
			//		Open the file using the usual search rules.
			// --------------------------------------------------------------------
			FileStream fp=Proj.pj_open_lib(ctx, gridname, FileAccess.Read);
			if(fp==null)
			{
				Proj.pj_errno=Libc.errno;
				ctx.last_errno=0; // don't treat as a persistent error
				return gilist;
			}

			gilist.filename=gridname;

			// --------------------------------------------------------------------
			//		Load a header, to determine the file type.
			// --------------------------------------------------------------------
			try
			{
				byte[] header=new byte[160];

				if(fp.Read(header, 0, header.Length)!=header.Length)
				{
					// some files may be smaller that sizeof(header), eg 160, so
					ctx.last_errno = 0; // don't treat as a persistent error
				}

				fp.Seek(0, SeekOrigin.Begin);

				// --------------------------------------------------------------------
				//		Determine file type.
				// --------------------------------------------------------------------
				if(Encoding.ASCII.GetString(header, 0, 6)=="HEADER"&&
					Encoding.ASCII.GetString(header, 96, 6)=="W GRID"&&
					Encoding.ASCII.GetString(header, 144, 16)=="TO      NAD83   ")
				{
					pj_gridinfo_init_ntv1(ctx, fp, gilist);
				}
				else if(Encoding.ASCII.GetString(header, 0, 8)=="NUM_OREC"&&
					Encoding.ASCII.GetString(header, 48, 7)=="GS_TYPE")
				{
					pj_gridinfo_init_ntv2(ctx, fp, gilist);
				}
				else if(gridname.Length>4&&gridname.EndsWith("gtx", StringComparison.CurrentCultureIgnoreCase))
				{
					pj_gridinfo_init_gtx(ctx, fp, gilist);
				}
				else if(Encoding.ASCII.GetString(header, 0, 9)=="CTABLE V2")
				{
					CTABLE ct=nad_ctable2_init(ctx, fp);

					gilist.format="ctable2";
					gilist.ct=ct;

					Proj.pj_log(ctx, PJ_LOG.DEBUG_MAJOR, "Ctable2 {0} {1}x{2}: LL=({3},{4}) UR=({5},{6})",
						ct.id, ct.lim.lam, ct.lim.phi, ct.ll.lam*Proj.RAD_TO_DEG, ct.ll.phi*Proj.RAD_TO_DEG,
						(ct.ll.lam+(ct.lim.lam-1)*ct.del.lam)*Proj.RAD_TO_DEG, (ct.ll.phi+(ct.lim.phi-1)*ct.del.phi)*Proj.RAD_TO_DEG);
				}
				else
				{
					CTABLE ct=nad_ctable_init(ctx, fp);

					if(ct==null)
					{
						Proj.pj_log(ctx, PJ_LOG.DEBUG_MAJOR, "CTABLE ct is NULL.");
					}
					else
					{
						gilist.format="ctable";
						gilist.ct=ct;

						Proj.pj_log(ctx, PJ_LOG.DEBUG_MAJOR, "Ctable {0} {1}x{2}: LL=({3},{4}) UR=({5},{6})",
							ct.id, ct.lim.lam, ct.lim.phi, ct.ll.lam*Proj.RAD_TO_DEG, ct.ll.phi*Proj.RAD_TO_DEG,
							(ct.ll.lam+(ct.lim.lam-1)*ct.del.lam)*Proj.RAD_TO_DEG, (ct.ll.phi+(ct.lim.phi-1)*ct.del.phi)*Proj.RAD_TO_DEG);
					}
				}
			}
			catch
			{
				Proj.pj_ctx_set_errno(ctx, -38);
			}

			fp.Close();

			return gilist;
		}
	}
}
