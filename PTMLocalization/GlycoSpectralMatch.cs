using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer.GlycoSearch
{
    public class GlycoSpectralMatch
    {
        public GlycoSpectralMatch(List<LocalizationGraph> localizationGraphs, GlycoType glycoType)
        {
            LocalizationGraphs = localizationGraphs;
            GlycanType = glycoType;
        }
        // placeholder
        public GlycoSpectralMatch()
        {

        }

        #region Proterties

        //Glyco properties
        public List<int> ModPos { get; set; }
        public LocalizationLevel LocalizationLevel { get; set; }

        public string localizerOutput { get; set; }

        public static GlycanBox GetFirstGraphGlycanBox(GlycoSpectralMatch gsm)
        {

            if (gsm.GlycanType == GlycoType.OGlycoPep)
            {
                return GlycanBox.OGlycanBoxes[gsm.LocalizationGraphs.First().ModBoxId];
            }
            else if (gsm.GlycanType == GlycoType.NGlycoPep)
            {
                return GlycanBox.NGlycanBoxes[gsm.LocalizationGraphs.First().ModBoxId];
            }
            else
            {
                return GlycanBox.MixedModBoxes[gsm.LocalizationGraphs.First().ModBoxId];
            }

        }
        
        public static GlycanBox GetGraphGlycanBox(LocalizationGraph graph, GlycoType type)
        {

            if (type == GlycoType.OGlycoPep)
            {
                return GlycanBox.OGlycanBoxes[graph.ModBoxId];
            }
            else if (type == GlycoType.NGlycoPep)
            {
                return GlycanBox.NGlycanBoxes[graph.ModBoxId];
            }
            else
            {
                return GlycanBox.MixedModBoxes[graph.ModBoxId];
            }

        }

        public static Glycan[] GetFirstGraphGlycans(GlycoSpectralMatch gsm, GlycanBox glycanBox)
        {
            var glycans = new Glycan[glycanBox.ModCount];
            if (gsm.GlycanType == GlycoType.OGlycoPep)
            {
                for (int i = 0; i < glycanBox.ModCount; i++)
                {
                    glycans[i] = GlycanBox.GlobalOGlycans[glycanBox.ModIds[i]];
                }
                return glycans;
            }
            else if (gsm.GlycanType == GlycoType.NGlycoPep)
            {
                for (int i = 0; i < glycanBox.ModCount; i++)
                {
                    glycans[i] = GlycanBox.GlobalNGlycans[glycanBox.ModIds[i]];
                }
                return glycans;
            }
            else
            {
                for (int i = 0; i < glycanBox.ModCount; i++)
                {
                    glycans[i] = GlycanBox.GlobalMixedGlycans[glycanBox.ModIds[i]];
                }
                return glycans;
            }
        }

        //Glycan type indicator
        public GlycoType GlycanType { get; } 

        public bool NGlycanMotifExist { get; set; } //NGlycan Motif exist. 


        //Glyco Info
        public List<LocalizationGraph> LocalizationGraphs { get; }  //Graph-based Localization information.
        public List<Route> Routes { get; set; } //Localized modification sites and modfication ID.
        public double DeltaScoreRoute { get; set; }  // difference between best and second-best routes in the best graph (i.e., localization delta score)
        public double DeltaScoreGraph { get; set; }  // difference between best routes in the best and second-best graphs (i.e., glycan composition delta score)
        public GlycanBox secondBestBox { get; set; }    // for recording next-best composition
        public double ScanInfo_p { get; set; }  //Scan P value, Used for Localization probability calculation. Ref PhosphoRS paper.

        public int Thero_n { get; set; } //Scan n value. Used for Localization probability calculation. Ref PhosphoRS paper.

        public Dictionary<int, List<Tuple<int, double>>> SiteSpeciLocalProb { get; set; } // Data <modPos, List<glycanId, site probability>>

        public List<GlycoSite> LocalizedGlycan { get; set; } 

        public double oxoRatio { get; set; }

        #endregion


        #region Glycopeptide Localization Output

        public static string AllLocalizationInfo(List<Route> routes)
        {
            string local = "";

            if (routes == null || routes.Count == 0)
            {
                return local;
            }
            //Some GSP have a lot paths, in which case only output first 10 paths and the total number of the paths.
            int maxOutputPath = 10;
            if (routes.Count <= maxOutputPath)
            {
                maxOutputPath = routes.Count;
            }

            int i = 0;
            while (i < maxOutputPath)
            {
                var ogl = routes[i];
                local += "{@" + ogl.ModBoxId.ToString() + "[";
                var g = string.Join(",", ogl.Mods.Select(p => (p.ModSite - 1).ToString() + "-" + p.GlycanID.ToString()));
                local += g + "]}";
                i++;
            }

            if (routes.Count > maxOutputPath)
            {
                local += "... In Total:" + routes.Count.ToString() + " Paths";
            }

            return local;
        }

        public static void LocalizedSiteSpeciLocalInfo(Dictionary<int, List<Tuple<int, double>>> siteSpeciLocalProb, List<GlycoSite> localizedGlycan, int? OneBasedStartResidueInProtein, Glycan[] globalGlycans, ref string local, ref string local_protein)
        {
            if (siteSpeciLocalProb == null)
            {
                return;
            }

            foreach (var loc in localizedGlycan.Where(p => p.IsLocalized))
            {
                var x = siteSpeciLocalProb[loc.ModSite].Where(p => p.Item1 == loc.GlycanID).First().Item2;
                var peptide_site = loc.ModSite - 1;
                local += "[" + peptide_site + "," + globalGlycans[loc.GlycanID].Composition + "," + x.ToString("0.000") + "]";

                //var protein_site = OneBasedStartResidueInProtein.HasValue ? OneBasedStartResidueInProtein.Value + loc.ModSite - 2 : -1;
                //local_protein += "[" + protein_site + "," + globalGlycans[loc.GlycanID].Composition + "," + x.ToString("0.000") + "]";
            }

        }
        
        public static string SiteSpeciLocalInfo(Dictionary<int, List<Tuple<int, double>>> siteSpeciLocalProb)
        {
            string local = "";

            if (siteSpeciLocalProb == null)
            {
                return local;
            }

            foreach (var sitep in siteSpeciLocalProb)
            {
                var site_1 = sitep.Key - 1;
                local += "{@" + site_1;
                foreach (var s in sitep.Value)
                {
                    local += "[" + s.Item1 + "," + s.Item2.ToString("0.000") + "]";
                }
                local += "}";
            }

            return local;
        }

        public Dictionary<int, List<GlycoSite>> GetGlycoSitesByGlycanID()
        {
            Dictionary<int, List<GlycoSite>> idDict = new();
            foreach (GlycoSite site in LocalizedGlycan)
            {
                if (idDict.ContainsKey(site.GlycanID))
                {
                    idDict[site.GlycanID].Add(site);
                } 
                else
                {
                    idDict[site.GlycanID] = new List<GlycoSite> {site};
                }
            }
            return idDict;
        }

        // calculate the difference in score between the top-scoring route of the best and second best graphs
        public static double CalculateGraphDeltaScore(GlycoSpectralMatch gsm, List<LocalizationGraph> allGraphs)
        {
            if (allGraphs.Count > 1)
            {
                var sortedGraphs = allGraphs.OrderByDescending(graph => graph.TotalScore).ToList();
                gsm.secondBestBox = GetGraphGlycanBox(sortedGraphs[1], gsm.GlycanType);
                return sortedGraphs[0].TotalScore - sortedGraphs[1].TotalScore;
            }

            // not applicable: need at least 2 graphs to calculate delta
            return -1;
        }

        #endregion

        public byte[] getTotalKind()
        {
            var glycanBox = GetFirstGraphGlycanBox(this);
            return glycanBox.Kind;
        }

        public string WriteLine(int? OneBasedStartResidueInProtein)
        {
            var sb = new StringBuilder();
            if (LocalizationGraphs != null)
            {
                sb.Append(LocalizationGraphs.First().TotalScore + "\t");

                var glycanBox = GetFirstGraphGlycanBox(this);

                sb.Append(glycanBox.ModCount + "\t");

                //sb.Append(LocalizationGraphs.First().ModPos.Length + "\t");

                //sb.Append(glycanBox.Mass + "\t");

                sb.Append(Glycan.GetFragpipeGlycan(glycanBox.Kind)); sb.Append("\t");

                //Get glycans
                var glycans = GetFirstGraphGlycans(this, glycanBox);

                sb.Append(string.Join(",", glycans.Select(p => p.Composition).ToArray()));
                sb.Append("\t");

                if (DeltaScoreGraph > -1)
                {
                    // leave blank if <2 graphs
                    sb.Append(DeltaScoreGraph);    
                    sb.Append('\t');
                    // get second graph compositions list
                    sb.Append(string.Join(",", GetFirstGraphGlycans(this, secondBestBox).Select(p => p.Composition).ToArray()));
                }
                else
                {
                    sb.Append('\t');
                }
                sb.Append('\t');

                if (Routes != null)
                {
                    sb.Append(LocalizationLevel); sb.Append("\t");
                    sb.Append(DeltaScoreRoute); sb.Append("\t");

                    string local_peptide = "";
                    string local_protein = "";
                    if (GlycanType == GlycoType.OGlycoPep)
                    {
                        LocalizedSiteSpeciLocalInfo(SiteSpeciLocalProb, LocalizedGlycan, OneBasedStartResidueInProtein, GlycanBox.GlobalOGlycans, ref local_peptide, ref local_protein);
                    }
                    else if (GlycanType == GlycoType.NGlycoPep)
                    {
                        LocalizedSiteSpeciLocalInfo(SiteSpeciLocalProb, LocalizedGlycan, OneBasedStartResidueInProtein, GlycanBox.GlobalNGlycans, ref local_peptide, ref local_protein);
                    }
                    else
                    {
                        LocalizedSiteSpeciLocalInfo(SiteSpeciLocalProb, LocalizedGlycan, OneBasedStartResidueInProtein, GlycanBox.GlobalMixedGlycans, ref local_peptide, ref local_protein);
                    }

                    sb.Append(local_peptide); sb.Append("\t");
                    //sb.Append(local_protein); sb.Append("\t");

                    //sb.Append(AllLocalizationInfo(Routes)); sb.Append("\t");

                    //sb.Append(SiteSpeciLocalInfo(SiteSpeciLocalProb)); sb.Append("\t");
                }
                else
                {
                    sb.Append(LocalizationLevel); sb.Append("\t");
                    //sb.Append("\t");
                    //sb.Append("\t");
                    //sb.Append("\t");
                    //sb.Append("\t");
                }
                sb.Append(oxoRatio); sb.Append("\t");
                sb.Append(NGlycanMotifExist);
            }
            return sb.ToString();
        }
    }
}