/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.bioalg.reportprocessor;

import java.awt.Point;
import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import javax.swing.JOptionPane;
import org.bioalg.graphprocessor.*;
import org.bioalg.imageprocessor.AImage;
import org.bioalg.imageprocessor.ImageIO;
import org.bioalg.mathprocessor.BasicMath;
import org.biolag.particleprocessor.Particle;
import org.biolag.particleprocessor.ParticleProcessor;
import org.biolag.particleprocessor.ParticlesModel;
import org.netbeans.modules.riceoptions.RiceOptionsOptionsPanelController;
import org.openide.util.NbPreferences;
import org.w3c.dom.Element;
import org.bioalg.imageprocessor.EDistanceMap;
import org.bioalg.imageprocessor.Kernels;

/**
 *
 * @author faroq
 */
public class Report {
    // main models

    private GraphModel gm;
    private ParticlesModel pm;
    private File imagefile;
    double numberofPrimaryBranches;
    double primaryBranches_lengthAverage;
    double secondaryBranches_lengthAverage;
    double numberofSecondaryBranches;
    double numberofTeratiaryBranches;
    double primary_int_average;
    double secondary_int_average;
    double rachis_length = 0.0; // the distance form the start to the end vertices
    double spikelet_nb = 0.0; // number of grains
    int numberofnodes = 0;
    ArrayList<BranchModel> primary_branches; // primary vertice positions
    AImage image;
    AImage subimage;
    AImage traimage;

    public Report(GraphModel gm) {
        this.gm = gm;
    }

    public Report() {
    }

    public Report(ParticlesModel pm) {
        this.pm = pm;
    }

    public Report(GraphModel gm, ParticlesModel pm) {
        this.gm = gm;
        this.pm = pm;

    }

    /**
     * Constructor by using xml files for the graph and particles models
     *
     * @param gmXML xml element of the graph stored in an xml file
     * @param pmXML xml element of the particles model stored in an xml file
     */
    public Report(Element gmXML, Element pmXML) {
        this.gm = new GraphModel(gmXML);
        this.pm = new ParticlesModel(pmXML);

    }

    public File getImagefile() {
        return imagefile;
    }

    public void setImagefile(File imagefile) {
        this.imagefile = imagefile;
    }

    /**
     * Collects the results and generate a master report
     *
     * @return
     */
    public String generateReport() {
        String report = "";
        if (gm != null) {
            // go for the graph model analyise it and collect results
            String mainstring = "";
            this.rachis_length = BasicMath.pixels2Inches(gm.getfixedPoints().get(0).distance(gm.getfixedPoints().get(1)));
            this.primary_branches = gm.getAllBranches();
            String branches = "";
            String branchesPositions = "";
            String displacements = "";
            String secondaryBranches = "";
            String secondaryBranchesPositions = "";
            String secondaryDisplacements = "";
            String tertiaryBranches = "";
            String tertiaryBranchesPositions = "";
            String tertiaryDisplacements = "";
            numberofPrimaryBranches = 0;
            numberofSecondaryBranches = 0;
            int count = 0;

            //first branch
            BranchModel bm = primary_branches.get(0);
            bm.calcpath();
            ArrayList<Path> paths = bm.getStart2endPath();

            for (int j = 0; j < paths.size(); j++) {
                Path path = paths.get(j);
                primaryBranches_lengthAverage += BasicMath.pixels2Inches(path.getlength());
                branches += "length_pb_" + 1 + "-" + (j + 1) + "," + BasicMath.pixels2Inches(path.getlength()) + "\n";
                branchesPositions += "position_pb" + 1 + "-" + (j + 1) + ",(x=" + BasicMath.pixels2Inches(path.get(0).x) + ": y=" + BasicMath.pixels2Inches(path.get(0).y) + ")\n";
                numberofPrimaryBranches++;
                count++;
                ArrayList<Path> secondaryPaths = bm.getSecondaryBranches(path);
                //secondary
                for (int k = 0; k < secondaryPaths.size(); k++) {
                    Path subpath = secondaryPaths.get(k);
                    secondaryBranches += "length_pb_" + 1 + "-" + (j + 1) + "_Sb" + (k + 1) + "," + BasicMath.pixels2Inches(subpath.getlength()) + "\n";
                    secondaryBranchesPositions += "position_pb" + 1 + "-" + (j + 1) + "_Sb" + (k + 1) + ",(x=" + BasicMath.pixels2Inches(subpath.get(0).x) + ": y=" + BasicMath.pixels2Inches(subpath.get(0).y) + ")\n";
                    numberofSecondaryBranches++;
                    secondaryDisplacements += "In_pb" + 1 + "-" + (j + 1) + "_Sb" + (k + 1) + "," + BasicMath.pixels2Inches(path.get(0).distance(subpath.get(0))) + "\n";
                    ArrayList<Path> tertiaryPaths = bm.getTertiaryBranches(subpath);
                    for (int l = 0; k < tertiaryPaths.size(); l++) {
                        ;
                    }
                }
            }


            // the rest of branches
            Vertex v1;
            if (paths.size() > 0 && paths.get(0).size() > 0) {
                v1 = paths.get(0).get(0);
            } else {
                return "paths.size():" + paths.size() + " Error: Unexpected error. Please check there are two main generating points in the structure model.\n";
            }

            Vertex v2 = null;
            for (int i = 1; i < primary_branches.size(); i++) {
                bm = primary_branches.get(i);
                bm.calcpath();
                paths = bm.getStart2endPath();

                for (int j = 0; j < paths.size(); j++) {
                    Path path = paths.get(j);
                    primaryBranches_lengthAverage += BasicMath.pixels2Inches(path.getlength());
                    branches += "length_pb_" + (i + 1) + "-" + (j + 1) + "," + BasicMath.pixels2Inches(path.getlength()) + "\n";
                    branchesPositions += "position_pb" + (i + 1) + "-" + (j + 1) + ",(x=" + BasicMath.pixels2Inches(path.get(0).x) + ": y=" + BasicMath.pixels2Inches(path.get(0).y) + ")\n";
                    numberofPrimaryBranches++;
                    count++;
                    v2 = paths.get(j).get(0);
                    ArrayList<Path> secondaryPaths = bm.getSecondaryBranches(path);
                    // secondary
                    for (int k = 0; k < secondaryPaths.size(); k++) {
                        Path subpath = secondaryPaths.get(k);
                        secondaryBranches += "length_pb_" + (i + 1) + "-" + (j + 1) + "_Sb" + (k + 1) + "," + BasicMath.pixels2Inches(subpath.getlength()) + "\n";
                        secondaryBranchesPositions += "position_pb" + (i + 1) + "-" + (j + 1) + "_Sb" + (k + 1) + ",(x=" + BasicMath.pixels2Inches(subpath.get(0).x) + ": y=" + BasicMath.pixels2Inches(subpath.get(0).y) + ")\n";
                        secondaryDisplacements += "In_pb" + (i + 1) + "-" + (j + 1) + "_Sb" + (k + 1) + "," + BasicMath.pixels2Inches(path.get(0).distance(subpath.get(0))) + "\n";
                        numberofSecondaryBranches++;
                    }

                }
                if (v2 != null) {
                    displacements += "In_pb" + (i) + "-pb" + (i + 1) + "," + BasicMath.pixels2Inches(v1.distance(v2)) + "\n";
                    v1 = v2;
                }

            }
            mainstring += "rachis_length, " + rachis_length + "\n"
                    + "pb_nb, " + numberofPrimaryBranches + "\n"
                    + "pb_Length_Average, " + BasicMath.roundTwoDecimals(primaryBranches_lengthAverage / numberofPrimaryBranches) + "\n"
                    + "pb_nb/rachis_lentgh," + BasicMath.roundTwoDecimals((double) (numberofPrimaryBranches) / (double) (rachis_length)) + "\n"
                    + "Primary Branches\n"
                    + branchesPositions + "\n"
                    + branches + "\n"
                    + displacements + "\n"
                    + "Secondary Branches\n"
                    + secondaryBranchesPositions + "\n"
                    + secondaryBranches + "\n"
                    + "Sb_nb," + numberofSecondaryBranches + "\n"
                    + "Sb_nb/Pb_nb," + BasicMath.roundTwoDecimals(numberofSecondaryBranches / numberofPrimaryBranches) + "\n"
                    + secondaryDisplacements;
            report = mainstring;
        }


        // THE PARTICLES
        if (pm != null) {
            //analyise the particles model
        }
        //System.out.println(report);
        return report;
    }

    public String generateReport3(String prevrecord) {


        String report = "";
        if (gm != null) {

            boolean calcgrain = (pm == null ? false : seperateJunctions());
            // go for the graph model analyise it and collect results
            String mainstring = "";
            // this.rachis_length = BasicMath.pixels2Inches(gm.getfixedPoints().get(0).distance(gm.getfixedPoints().get(1)));
            this.primary_branches = gm.getAllBranches();
            String branches = "";
            String branchesPositions = "";
            String displacements = "";
            String secondaryBranches = "";
            String secondaryBranchesPositions = "";
            String secondaryDisplacements = "";
            String tertiaryBranches = "";
            String tertiaryBranchesPositions = "";
            String tertiaryDisplacements = "";
            numberofPrimaryBranches = 0;
            numberofSecondaryBranches = 0;
            numberofnodes = 1;
            String record = "";
            int count = 0;

            if (primary_branches.isEmpty()) {
                return "\n";
            }
            //first branch
            BranchModel bm = primary_branches.get(0);
            bm.calcpath();
            ArrayList<Path> paths = bm.getStart2endPath();

            for (int j = 0; j < paths.size(); j++) {
                Path path = paths.get(j);

                // pb length node
                branches = "SA" + (count + 1) + "," + BasicMath.pixels2Inches(path.getlength()) + "," + numberofnodes;


                count++;
                ArrayList<Path> secondaryPaths = bm.getSecondaryBranches(path);
                //secondary

                // sb_nb sp_np interval
                branches += "," + secondaryPaths.size() + ","
                        + // grains
                        (calcgrain == true ? ParticleProcessor.ispariclesinArea(image, path.get(path.size() - 1), path.get(0), pm.getparticles()) : 0)
                        +",0,";
                for (int k = 0; k < secondaryPaths.size(); k++) {
                    Path subpath = secondaryPaths.get(k);
                    ArrayList<Path> tertiaryPaths = bm.getTertiaryBranches(subpath);
                    // sp length node tb_nb sp_nb interval
                    secondaryBranches = "SA" + (count) + "_TA" + (k + 1) + "," + BasicMath.pixels2Inches(subpath.getlength())
                             + "," + tertiaryPaths.size() + ","
                            + (calcgrain == true ? ParticleProcessor.ispariclesinArea(subimage, subpath.get(subpath.size() - 1), subpath.get(0), pm.getparticles()) : 0)
                            + "," + BasicMath.pixels2Inches(path.get(0).distance(subpath.get(0)));
                            

                    secondaryBranches += ",";
                    for (int l = 0; l < tertiaryPaths.size(); l++) {
                        Path terpath = tertiaryPaths.get(l);
                        tertiaryBranches = "SA" + (count) + "_TA" + (k + 1) + "_Qb" + (l + 1)
                                + "," + BasicMath.pixels2Inches(terpath.getlength()) + ","
                                + (calcgrain == true ? ParticleProcessor.ispariclesinArea(traimage, terpath.get(terpath.size() - 1), terpath.get(0), pm.getparticles()) : 0)
                                + "," + BasicMath.pixels2Inches(subpath.get(0).distance(terpath.get(0)));
                        record += prevrecord + branches + secondaryBranches + tertiaryBranches + "\n";
                    }
                    if (tertiaryPaths.size() <= 0) {
                        record += prevrecord + branches + secondaryBranches + "NV,NV,NV,NV\n";
                    }

                }
                if (secondaryPaths.isEmpty()) {
                    record += prevrecord + branches + "NV,NV,NV,NV,NV,NV,NV,NV,NV,NV\n";
                }
            }


            // the rest of branches
            Vertex v1;
            if (paths.size() > 0 && paths.get(0).size() > 0) {
                v1 = paths.get(0).get(0);
            } else {
                return "paths.size():" + paths.size() + " Error: Please check start and end cycles are existing and/or no cycles are in the structure.\n";
            }


            Vertex v2 = null;
            for (int i = 1; i < primary_branches.size(); i++) {
                bm = primary_branches.get(i);
                bm.calcpath();
                paths = bm.getStart2endPath();

                for (int j = 0; j < paths.size(); j++) {
                    Path path = paths.get(j);
                    v2 = paths.get(j).get(0);
                    double dist = 0; // for interval
                    if (v2 != null) {
                        String minnodesDistance = NbPreferences.forModule(RiceOptionsOptionsPanelController.class).get("nodemindistance", "");
                        double minnodedist = (BasicMath.isdouble(minnodesDistance) ? Double.parseDouble(minnodesDistance) : 10);
                        dist = BasicMath.pixels2Inches(v1.distance(v2));
                        if (v1.distance(v2) > minnodedist) {

                            numberofnodes++;

                        }
                        v1 = v2;
                    }
                    // pb length node
                    branches = "SA" + (count + 1) + "," + BasicMath.pixels2Inches(path.getlength()) + "," + (numberofnodes);
                    count++;

                    ArrayList<Path> secondaryPaths = bm.getSecondaryBranches(path);
                    branches += "," + secondaryPaths.size() + ","
                            + (calcgrain == true ? ParticleProcessor.ispariclesinArea(image, path.get(path.size() - 1), path.get(0), pm.getparticles()) : 0)
                            + "," + dist + ",";


                    // secondary
                    for (int k = 0; k < secondaryPaths.size(); k++) {
                        Path subpath = secondaryPaths.get(k);
                        ArrayList<Path> tertiaryPaths = bm.getTertiaryBranches(subpath);
                        // sp length node tb_nb sp_nb interval
                        secondaryBranches = "SA" + (count) + "_TA" + (k + 1) + "," + BasicMath.pixels2Inches(subpath.getlength())
                                + "," + tertiaryPaths.size() + ","
                                + (calcgrain == true ? ParticleProcessor.ispariclesinArea(subimage, subpath.get(subpath.size() - 1), subpath.get(0), pm.getparticles()) : 0)
                                + "," + BasicMath.pixels2Inches(path.get(0).distance(subpath.get(0)));
                        secondaryBranches += ",";
                        for (int l = 0; l < tertiaryPaths.size(); l++) {
                            Path terpath = tertiaryPaths.get(l);
                            tertiaryBranches = "SA" + (count) + "_TA" + (k + 1) + "_QA" + (l + 1)
                                    + "," + BasicMath.pixels2Inches(terpath.getlength()) + ","
                                    + (calcgrain == true ? ParticleProcessor.ispariclesinArea(traimage, terpath.get(terpath.size() - 1), terpath.get(0), pm.getparticles()) : 0)
                                    + "," + BasicMath.pixels2Inches(subpath.get(0).distance(terpath.get(0)));
                            record += prevrecord + branches + secondaryBranches + tertiaryBranches + "\n";
                        }
                        if (tertiaryPaths.size() <= 0) {
                            record += prevrecord + branches + secondaryBranches + "NV,NV,NV,NV\n";
                        }
                    }
                    if (secondaryPaths.isEmpty()) {
                        record += prevrecord + branches + "NV,NV,NV,NV,NV,NV,NV,NV,NV,NV\n";
                    }

                }


            }
            mainstring += record;
            report = mainstring;
        }


        // THE PARTICLES
        if (pm != null) {
            //analyise the particles model
        }
        //System.out.println(report);
        return report;
    }

    public GraphModel getGraphModel() {
        return gm;
    }

    public void setGm(GraphModel gm) {
        this.gm = gm;
    }

    public ParticlesModel getPm() {
        return pm;
    }

    public void setPm(ParticlesModel pm) {
        this.pm = pm;
    }

    public String generateReport2() {


        String report = "";
        if (gm != null) {
            // go for the graph model analyise it and collect results
            String mainstring = "";
            //this.rachis_length = BasicMath.pixels2Inches(gm.getfixedPoints().get(0).distance(gm.getfixedPoints().get(1)));
            this.primary_branches = gm.getAllBranches();
            this.rachis_length = BasicMath.pixels2Inches(gm.getRachisLength(primary_branches));
            numberofnodes = 0;
            String branches = "";
            String branchesPositions = "";
            String displacements = "";
            String secondaryBranches = "";
            String secondaryBranchesPositions = "";
            String secondaryDisplacements = "";
            double pb_int_average = 0;
            double sb_int_average = 0;
            int count = 0;
            numberofPrimaryBranches = 0;
            numberofSecondaryBranches = 0;
            //first branch
            if (!primary_branches.isEmpty()) {
                BranchModel bm = primary_branches.get(0);
                bm.calcpath();
                ArrayList<Path> paths = bm.getStart2endPath();

                for (int j = 0; j < paths.size(); j++) {
                    Path path = paths.get(j);
                    primaryBranches_lengthAverage += BasicMath.pixels2Inches(path.getlength());
                    branches += "length_pb_" + 1 + "-" + (j + 1) + "," + BasicMath.pixels2Inches(path.getlength()) + "\n";
                    branchesPositions += "position_pb" + 1 + "-" + (j + 1) + ",(x=" + BasicMath.pixels2Inches(path.get(0).x) + ": y=" + BasicMath.pixels2Inches(path.get(0).y) + ")\n";
                    numberofPrimaryBranches++;
                    count++;
                    ArrayList<Path> secondaryPaths = bm.getSecondaryBranches(path);
                    //secondary
                    for (int k = 0; k < secondaryPaths.size(); k++) {
                        Path subpath = secondaryPaths.get(k);
                        secondaryBranches += "length_pb_" + 1 + "-" + (j + 1) + "_Sb" + (k + 1) + "," + BasicMath.pixels2Inches(subpath.getlength()) + "\n";
                        secondaryBranchesPositions += "position_pb" + 1 + "-" + (j + 1) + "_Sb" + (k + 1) + ",(x=" + BasicMath.pixels2Inches(subpath.get(0).x) + ": y=" + BasicMath.pixels2Inches(subpath.get(0).y) + ")\n";
                        numberofSecondaryBranches++;
                        secondaryDisplacements += "ln_pb" + 1 + "-" + (j + 1) + "_Sb" + (k + 1) + "," + BasicMath.pixels2Inches(path.get(0).distance(subpath.get(0))) + "\n";
                        sb_int_average += BasicMath.pixels2Inches(path.get(0).distance(subpath.get(0)));

                        ArrayList<Path> tertiaryPaths = bm.getTertiaryBranches(subpath);
                        for (int l = 0; l < tertiaryPaths.size(); l++) {
                            numberofTeratiaryBranches++;
                        }
                    }
                }


                // the rest of branches

                numberofnodes = 1;
                Vertex v1;
                if (paths.size() > 0 && paths.get(0).size() > 0) {
                    v1 = paths.get(0).get(0);
                } else {
                    return "paths.size():" + paths.size() + " Error: Please check start and end cycles are exists and/or no cycles are in the structure.\n";
                }
                Vertex v2 = null;
                for (int i = 1; i < primary_branches.size(); i++) {
                    bm = primary_branches.get(i);
                    bm.calcpath();
                    paths = bm.getStart2endPath();

                    for (int j = 0; j < paths.size(); j++) {
                        Path path = paths.get(j);
                        primaryBranches_lengthAverage += BasicMath.pixels2Inches(path.getlength());
                        branches += "length_pb_" + (i + 1) + "-" + (j + 1) + "," + BasicMath.pixels2Inches(path.getlength()) + "\n";
                        branchesPositions += "position_pb" + (i + 1) + "-" + (j + 1) + ",(x=" + BasicMath.pixels2Inches(path.get(0).x) + ": y=" + BasicMath.pixels2Inches(path.get(0).y) + ")\n";
                        numberofPrimaryBranches++;
                        count++;
                        v2 = paths.get(j).get(0);
                        ArrayList<Path> secondaryPaths = bm.getSecondaryBranches(path);
                        // secondary
                        for (int k = 0; k < secondaryPaths.size(); k++) {
                            Path subpath = secondaryPaths.get(k);
                            secondaryBranches += "length_pb_" + (i + 1) + "-" + (j + 1) + "_Sb" + (k + 1) + "," + BasicMath.pixels2Inches(subpath.getlength()) + "\n";
                            secondaryBranchesPositions += "position_pb" + (i + 1) + "-" + (j + 1) + "_Sb" + (k + 1) + ",(x=" + BasicMath.pixels2Inches(subpath.get(0).x) + ": y=" + BasicMath.pixels2Inches(subpath.get(0).y) + ")\n";
                            secondaryDisplacements += "ln_pb" + (i + 1) + "-" + (j + 1) + "_Sb" + (k + 1) + "," + BasicMath.pixels2Inches(path.get(0).distance(subpath.get(0))) + "\n";
                            sb_int_average += BasicMath.pixels2Inches(path.get(0).distance(subpath.get(0)));
                            numberofSecondaryBranches++;
                            ArrayList<Path> tertiaryPaths = bm.getTertiaryBranches(subpath);
                            for (int l = 0; l < tertiaryPaths.size(); l++) {
                                numberofTeratiaryBranches++;
                            }
                        }

                    }
                    if (v2 != null) {
                        displacements += "ln_pb" + (i) + "-pb" + (i + 1) + "," + BasicMath.pixels2Inches(v1.distance(v2)) + "\n";
                        pb_int_average += BasicMath.pixels2Inches(v1.distance(v2));

                        // read the minimum distance between the nodes
                        String minnodesDistance = NbPreferences.forModule(RiceOptionsOptionsPanelController.class).get("nodemindistance", "");
                        double minnodedist = (BasicMath.isdouble(minnodesDistance) ? Double.parseDouble(minnodesDistance) : 10);
                        if (v1.distance(v2) > minnodedist) {
                            numberofnodes++;

                        }
                        v1 = v2;
                    }

                }

                String execludenodes = NbPreferences.forModule(RiceOptionsOptionsPanelController.class).get("execludeNode", "");
                int execlude = (BasicMath.isInteger(execludenodes) ? Integer.parseInt(execludenodes) : 0);
                System.out.println(execludenodes);
                mainstring +=
                        rachis_length + ","
                        + 2 * BasicMath.pixels2Inches(getPanicleRadius(gm.getfixedPoints().get(0), 20)) + "," // 20 is the radius of the search                           
                        + (execlude != 0 ? numberofnodes - 2 : numberofnodes) + ","
                        + numberofPrimaryBranches + ","
                        + BasicMath.roundTwoDecimals(primaryBranches_lengthAverage / numberofPrimaryBranches) + ","
                        //  + BasicMath.roundTwoDecimals((double) (numberofPrimaryBranches) / (double) (rachis_length)) + ","
                        + BasicMath.roundTwoDecimals(pb_int_average / (numberofPrimaryBranches - 1)) + ","
                        + numberofSecondaryBranches + ","
                        + BasicMath.roundTwoDecimals(numberofSecondaryBranches / numberofPrimaryBranches) + ","
                        + BasicMath.roundTwoDecimals((sb_int_average / (numberofSecondaryBranches - 1))) + ","
                        + numberofTeratiaryBranches;
                report = mainstring;
            }
        }
        // FIND PARTICLES
        if (pm != null) {
            //analyise the particles model
            if (report.equals("")) {
                report += "0,0,0,0,0,0,0,0,0,0";
            }
            report += "," + String.valueOf(pm.getProportionalitySum());
//                    + "," + String.valueOf(BasicMath.pixels2Inches(pm.getwidthAverage()))
//                    + "," + String.valueOf(BasicMath.pixels2Inches(pm.getheightAverage()))
//                    + "," + String.valueOf(BasicMath.pixels2Inches(pm.getarea()))
//                    + "," + String.valueOf(BasicMath.pixels2Inches(pm.getprimeter()))
//                    + "," + String.valueOf(BasicMath.roundTwoDecimals(pm.getcircularity()))
//                    + "," + String.valueOf(BasicMath.roundTwoDecimals(pm.getcompactness()));

        }

        //System.out.println(report);
        return report;
    }

    public String GetParticlesReport() {
        String header = "ID,X,Y,No.Grains,Length,Width,Area,Primeter,Circularity,Compactness,Ellipticity,AR";
        String report = "";

        if (pm != null) {
            //pm.calcAll();
            ArrayList<Particle> particles = pm.getparticles();
            int size = particles.size();
            for (int i = 0; i < size; i++) {
                Particle particle = particles.get(i);
                report += i
                        + "," + particle.getPosition().x
                        + "," + particle.getPosition().y
                        + "," + String.valueOf(particle.getProportionality())
                        + "," + String.valueOf(BasicMath.pixels2Inches(particle.getWidth()))
                        + "," + String.valueOf(BasicMath.pixels2Inches(particle.getHeight()))
                        + "," + String.valueOf(BasicMath.pixels2Inches(particle.getArea() ))
                        + "," + String.valueOf(BasicMath.pixels2Inches(particle.getPrimeter() ))
                        + "," + String.valueOf(BasicMath.roundTwoDecimals(particle.getCircularity())) //String.valueOf(BasicMath.roundTwoDecimals(particle.getCircularity() > 1 ? 1 : particle.getCircularity()))
                        + "," + String.valueOf(BasicMath.roundTwoDecimals(particle.getCompactness()))
                        + "," + String.valueOf(BasicMath.roundTwoDecimals(particle.getEllipticity()))
                        + "," + String.valueOf(BasicMath.roundTwoDecimals(((double)particle.getWidth()) /((double) particle.getHeight())))
                        + "\n";

            }
        }
        return header + "\n" + report;


    }

    public String GetSumParticlesReport() {

        String report = "";

        if (pm != null) {
            //pm.calcAll();
            if (pm != null) {
                //analyise the particles model
                if (report.equals("")) {
                }
                report += String.valueOf(pm.getProportionalitySum())
                        + "," + String.valueOf(BasicMath.pixels2Inches(pm.getwidthAverage()))
                        + "," + String.valueOf(BasicMath.pixels2Inches(pm.getheightAverage()))
                        + "," + String.valueOf(BasicMath.pixels2Inches(pm.getarea()))
                        + "," + String.valueOf(BasicMath.pixels2Inches(pm.getprimeter()))
                        + "," + String.valueOf(BasicMath.roundTwoDecimals(pm.getcircularity()))
                        + "," + String.valueOf(BasicMath.roundTwoDecimals(pm.getcompactness()))
                        + "," + String.valueOf(BasicMath.roundTwoDecimals(pm.getellipciticiy()))
                        + "," + String.valueOf(BasicMath.roundTwoDecimals(pm.getAR()));


            }



        }
        return report;
    }

    //  helper function
    public boolean seperateJunctions() {
        if (imagefile == null) {
            JOptionPane.showMessageDialog(null, "Can not find the image file", "Report Generating", JOptionPane.ERROR_MESSAGE);
            return false;
        }

        image = ImageIO.getImageAsGrayscaledAImage(imagefile);
        if (image == null) {
            return false;
        }

        image.gaussianBlur(3, 1.0 / 9); // gaussian
        image = image.fastMeanThresholding(77, 5);
        subimage = new AImage(image);
        traimage = new AImage(image);

        int bigk = 25;
        int medk = 11;


        for (int i = 0; i < gm.getEdges().size(); i++) {
            EdgeModel edgeModel = gm.getEdges().get(i);
            Point p1 = new Point(edgeModel.getVertex1().getLocation());
            Point p2 = new Point(edgeModel.getVertex2().getLocation());

            p1 = BasicMath.swapXY(p1);
            p2 = BasicMath.swapXY(p2);

//            if (edgeModel.getVertex1().getType() == Vertex.VERTEX_TYPE.Tertiary) {
//                traimage.putMark(p1, medk, 0);
//            }
//            if (edgeModel.getVertex2().getType() == Vertex.VERTEX_TYPE.Tertiary) {
//                traimage.putMark(p2, medk, 0);
//            }
//            if (edgeModel.getVertex1().getType() == Vertex.VERTEX_TYPE.Seconday) {
//                subimage.putMark(p1, medk, 0);
//            }
//            if (edgeModel.getVertex2().getType() == Vertex.VERTEX_TYPE.Seconday) {
//                subimage.putMark(p2, medk, 0);
//            }

            if (edgeModel.getVertex1().getType() == Vertex.VERTEX_TYPE.Primary || edgeModel.getVertex1().getType() == Vertex.VERTEX_TYPE.Generating) {
                image.putMark(p1, bigk, 0);
            }
            if (edgeModel.getVertex2().getType() == Vertex.VERTEX_TYPE.Primary || edgeModel.getVertex2().getType() == Vertex.VERTEX_TYPE.Generating) {
                image.putMark(p2, bigk, 0);
            }
        }
        //ImageIO.saveAImageAsImage(image, "/home/faroq/comp.tiff");
        subimage = new AImage(image);
        traimage = new AImage(image);
        return true;
    }

    /**
     * find the panicle diameter by using the Euclidian Distance Map within a
     * search range centered at the first generating point (center parameter).
     * The diameter is just the largest radius of the maximal disk inside the
     * search
     *
     * @param center the point where the search has to be made
     * @param searchRadius the search radius centered at the center point
     * @return
     */
    public double getPanicleRadius(Point center, int searchRadius) {
        double panicleRadius = 0.0f;
        if (imagefile != null) {
            image = ImageIO.getImageAsGrayscaledAImage(imagefile);
            image = image.meanLowpass(5);
            image = image.fastMeanThresholding(75, 5);
        }
        // get the EDM image
        AImage EDMImage = EDistanceMap.meijster(image);
        //ImageIO.saveAImageAsImage(EDMImage, "/home/faroq/edm.tiff");

        // make a disk kernel of radius searchRadius center at the center point
        ArrayList<Point> searchSpace = Kernels.generateDiskAsPoints(searchRadius, center);
        // search and find the maximum 
        Iterator pointsIterator = searchSpace.iterator();
        while (pointsIterator.hasNext()) {
            Point currentPoint = (Point) pointsIterator.next();
            double currentPoint_radius = Math.sqrt(EDMImage.getRC(currentPoint.y, currentPoint.x, 0));

            if (currentPoint_radius > panicleRadius) {
                panicleRadius = currentPoint_radius;
            }
        }
        //JOptionPane.showMessageDialog(null,  "currentPoint_radius"+ panicleRadius + "\n Centerted at:"+ center.toString() );

        return panicleRadius;
    }
}
