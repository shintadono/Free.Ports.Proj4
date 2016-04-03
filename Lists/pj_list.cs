using Free.Ports.Proj4.Projections;

namespace Free.Ports.Proj4
{
	public static partial class Proj
	{
		public static PJ GetPJ(string name)
		{
			if(name==null||name.Length==0) return null;
			switch(name[0])
			{
				case 'a':
					switch(name)
					{
						case "aea": return new PJ_aea();
						case "aeqd": return new PJ_aeqd();
						case "airy": return new PJ_airy();
						case "aitoff": return new PJ_aitoff();
						case "alsk": return new PJ_alsk();
						case "apian": return new PJ_apian();
						case "august": return new PJ_august();
						default: return null;
					}
				case 'b':
					switch(name)
					{
						case "bacon": return new PJ_bacon();
						case "bipc": return new PJ_bipc();
						case "boggs": return new PJ_boggs();
						case "bonne": return new PJ_bonne();
						default: return null;
					}
				case 'c':
					switch(name)
					{
						case "calcofi": return new PJ_calcofi();
						case "cass": return new PJ_cass();
						case "cc": return new PJ_cc();
						case "cea": return new PJ_cea();
						case "chamb": return new PJ_chamb();
						case "collg": return new PJ_collg();
						case "crast": return new PJ_crast();
						default: return null;
					}
				case 'd':
					return name=="denoy"?new PJ_denoy():null;
				case 'e':
					switch(name)
					{
						case "eck1": return new PJ_eck1();
						case "eck2": return new PJ_eck2();
						case "eck3": return new PJ_eck3();
						case "eck4": return new PJ_eck4();
						case "eck5": return new PJ_eck5();
						case "eck6": return new PJ_eck6();
						case "eqc": return new PJ_eqc();
						case "eqdc": return new PJ_eqdc();
						case "etmerc": return new PJ_etmerc();
						case "euler": return new PJ_euler();
						default: return null;
					}
				case 'f':
					switch(name)
					{
						case "fahey": return new PJ_fahey();
						case "fouc": return new PJ_fouc();
						case "fouc_s": return new PJ_fouc_s();
						case "ftmerc": return new PJ_ftmerc();
						default: return null;
					}
				case 'g':
					switch(name)
					{
						case "gall": return new PJ_gall();
						case "geocent": return new PJ_geocent();
						case "geos": return new PJ_geos();
						case "gins8": return new PJ_gins8();
						case "gn_sinu": return new PJ_gn_sinu();
						case "gnom": return new PJ_gnom();
						case "goode": return new PJ_goode();
						case "gs48": return new PJ_gs48();
						case "gs50": return new PJ_gs50();
						case "gstmerc": return new PJ_gstmerc();
						default: return null;
					}
				case 'h':
					switch(name)
					{
						case "hammer": return new PJ_hammer();
						case "hatano": return new PJ_hatano();
						case "healpix": return new PJ_healpix();
						default: return null;
					}
				case 'i':
					switch(name)
					{
						case "igh": return new PJ_igh(); // "Interrupted Goode Homolosine"
						case "isea": return new PJ_isea();
						case "imw_p": return new PJ_imw_p();
						default: return null;
					}
				case 'k':
					switch(name)
					{
						case "kav5": return new PJ_kav5();
						case "kav7": return new PJ_kav7();
						case "krovak": return new PJ_krovak();
						case "ktmerc": return new PJ_ktmerc();
						default: return null;
					}
				case 'l':
					switch(name)
					{
						case "labrd": return new PJ_labrd();
						case "laea": return new PJ_laea();
						case "lagrng": return new PJ_lagrng();
						case "larr": return new PJ_larr();
						case "lask": return new PJ_lask();
						case "latlon": return new PJ_latlon();
						case "latlong": return new PJ_latlong();
						case "latlong_deg": return new PJ_latlong_deg(); // Add by the Authors of the C# Port
						case "lcc": return new PJ_lcc();
						case "lcca": return new PJ_lcca();
						case "leac": return new PJ_leac();
						case "lee_os": return new PJ_lee_os();
						case "longlat": return new PJ_longlat();
						case "lonlat": return new PJ_lonlat();
						case "loxim": return new PJ_loxim();
						case "lsat": return new PJ_lsat();
						default: return null;
					}
				case 'm':
					switch(name)
					{
						case "mbt_s": return new PJ_mbt_s();
						case "mbt_fps": return new PJ_mbt_fps();
						case "mbtfpp": return new PJ_mbtfpp();
						case "mbtfpq": return new PJ_mbtfpq();
						case "mbtfps": return new PJ_mbtfps();
						case "merc": return new PJ_merc();
						case "mil_os": return new PJ_mil_os();
						case "mill": return new PJ_mill();
						case "misrsom": return new PJ_misrsom();
						case "moll": return new PJ_moll();
						case "murd1": return new PJ_murd1();
						case "murd2": return new PJ_murd2();
						case "murd3": return new PJ_murd3();
						default: return null;
					}
				case 'n':
					switch(name)
					{
						case "natearth": return new PJ_natearth();
						case "nell": return new PJ_nell();
						case "nell_h": return new PJ_nell_h();
						case "nicol": return new PJ_nicol();
						case "nsper": return new PJ_nsper();
						case "nzmg": return new PJ_nzmg();
						default: return null;
					}
				case 'o':
					switch(name)
					{
						case "ob_tran": return new PJ_ob_tran();
						case "ocea": return new PJ_ocea();
						case "oea": return new PJ_oea();
						case "omerc": return new PJ_omerc();
						case "ortel": return new PJ_ortel();
						case "ortho": return new PJ_ortho();
						default: return null;
					}
				case 'p':
					switch(name)
					{
						case "patterson": return new PJ_patterson();
						case "pconic": return new PJ_pconic();
						case "poly": return new PJ_poly();
						case "putp1": return new PJ_putp1();
						case "putp2": return new PJ_putp2();
						case "putp3": return new PJ_putp3();
						case "putp3p": return new PJ_putp3p();
						case "putp4p": return new PJ_putp4p();
						case "putp5": return new PJ_putp5();
						case "putp5p": return new PJ_putp5p();
						case "putp6": return new PJ_putp6();
						case "putp6p": return new PJ_putp6p();
						default: return null;
					}
				case 'q':
					switch(name)
					{
						case "qua_aut": return new PJ_qua_aut();
						case "qsc": return new PJ_qsc();
						default: return null;
					}
				case 'r':
					switch(name)
					{
						case "rhealpix": return new PJ_rhealpix();
						case "robin": return new PJ_robin();
						case "rouss": return new PJ_rouss();
						case "rpoly": return new PJ_rpoly();
						default: return null;
					}
				case 's':
					switch(name)
					{
						case "sch": return new PJ_sch();
						case "sinu": return new PJ_sinu();
						case "somerc": return new PJ_somerc();
						case "stere": return new PJ_stere();
						case "sterea": return new PJ_sterea();
						default: return null;
					}
				case 't':
					switch(name)
					{
						case "tcc": return new PJ_tcc();
						case "tcea": return new PJ_tcea();
						case "tissot": return new PJ_tissot();
						case "tmerc": return new PJ_tmerc();
						case "tpeqd": return new PJ_tpeqd();
						case "tpers": return new PJ_tpers();
						default: return null;
					}
				case 'u':
					switch(name)
					{
						case "ups": return new PJ_ups();
						case "urm5": return new PJ_urm5();
						case "urmfps": return new PJ_urmfps();
						case "utm": return new PJ_utm();
						default: return null;
					}
				case 'v':
					switch(name)
					{
						case "vandg": return new PJ_vandg();
						case "vandg2": return new PJ_vandg2();
						case "vandg3": return new PJ_vandg3();
						case "vandg4": return new PJ_vandg4();
						case "vitk1": return new PJ_vitk1();
						default: return null;
					}
				case 'w':
					switch(name)
					{
						case "wag1": return new PJ_wag1();
						case "wag2": return new PJ_wag2();
						case "wag3": return new PJ_wag3();
						case "wag4": return new PJ_wag4();
						case "wag5": return new PJ_wag5();
						case "wag6": return new PJ_wag6();
						case "wag7": return new PJ_wag7();
						case "weren": return new PJ_weren();
						case "wink1": return new PJ_wink1();
						case "wink2": return new PJ_wink2();
						case "wintri": return new PJ_wintri();
						default: return null;
					}
				default: return null;
			}
		}
	}
}
