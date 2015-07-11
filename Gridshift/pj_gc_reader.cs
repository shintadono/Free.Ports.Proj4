//*****************************************************************************
//
// Project:  PROJ.4
// Purpose:  Code to read a grid catalog from a .cvs file.
// Author:   Frank Warmerdam, warmerdam@pobox.com
//
//*****************************************************************************
// Copyright (c) 2012, Frank Warmerdam <warmerdam@pobox.com>
// Copyright (c) 2014-2015 by the Authors
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

namespace Free.Ports.Proj4.Gridshift
{
	public static partial class Grid
	{
		//***********************************************************************
		//							pj_gc_readcatalog()
		//
		//		Read a grid catalog from a .csv file.
		//***********************************************************************
		public static PJ_GridCatalog pj_gc_readcatalog(projCtx ctx, string catalog_name)
		{
			FileStream fs=Proj.pj_open_lib(ctx, catalog_name, FileAccess.Read);

			if(fs==null) return null;

			using(StreamReader fid=new StreamReader(fs))
			{
				// discard title line
				fid.ReadLine();

				PJ_GridCatalog catalog;

				try
				{
					catalog=new PJ_GridCatalog();
				}
				catch
				{
					return null;
				}

				catalog.catalog_name=catalog_name;

				catalog.entries=new List<PJ_GridCatalog.Entry>();

				PJ_GridCatalog.Entry entry;
				while(!pj_gc_readentry(ctx, fid, out entry)) catalog.entries.Add(entry);

				fid.Close();

				return catalog;
			}
		}

		//***********************************************************************
		//							pj_gc_read_csv_line()
		//
		//		Simple csv line splitter with fixed maximum line size and
		//		token count.
		//***********************************************************************
		static char[] csv_line_splitter=new char[] { ',' };
		static string[] pj_gc_read_csv_line(projCtx ctx, StreamReader fid)
		{
			while(!fid.EndOfStream)
			{
				string line=fid.ReadLine();
				if(line==null) break;

				line=line.TrimStart();

				int earlyNULChar=line.IndexOf('\0');
				if(earlyNULChar>=0) line=line.Substring(0, earlyNULChar);

				// skip blank and comment lines
				if(line.Length==0) continue;
				if(line[0]=='#') continue;

				return line.Split(csv_line_splitter, StringSplitOptions.RemoveEmptyEntries);
			}

			return null;
		}

		//***********************************************************************
		//								pj_gc_parsedate()
		//
		//		Parse a date into a floating point year value. Acceptable
		//		values are "yyyy.fraction" and "yyyy-mm-dd". Anything else
		//		returns 0.0.
		//***********************************************************************
		public static double pj_gc_parsedate(projCtx ctx, string date_string)
		{
			return Proj.pj_parsedate(ctx, date_string);
		}

		//***********************************************************************
		//							pj_gc_readentry()
		//
		//		Read one catalog entry from the file.
		//
		//		Format:
		//			gridname,ll_long,ll_lat,ur_long,ur_lat,priority,date
		//***********************************************************************
		static bool pj_gc_readentry(projCtx ctx, StreamReader fid, out PJ_GridCatalog.Entry entry)
		{
			entry=new PJ_GridCatalog.Entry();

			string[] tokens=pj_gc_read_csv_line(ctx, fid);
			if(tokens==null||tokens.Length<5)
			{
				if(tokens.Length!=0) Proj.pj_log(ctx, PJ_LOG.ERROR, "Short line in grid catalog.");
				return true;
			}

			entry.definition=tokens[0];
			entry.region.ll_long=Proj.dmstor_ctx(ctx, tokens[1]);
			entry.region.ll_lat=Proj.dmstor_ctx(ctx, tokens[2]);
			entry.region.ur_long=Proj.dmstor_ctx(ctx, tokens[3]);
			entry.region.ur_lat=Proj.dmstor_ctx(ctx, tokens[4]);
			if(tokens.Length>5) int.TryParse(tokens[5], out entry.priority); // defaults to zero
			if(tokens.Length>6) entry.date=pj_gc_parsedate(ctx, tokens[6]);

			return false;
		}
	}
}
