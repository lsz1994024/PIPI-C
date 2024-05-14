package proteomics.PTM;

import ProteomicsLibrary.MassTool;
import gurobi.*;
import org.apache.commons.math3.util.Pair;
import org.dom4j.Document;
import org.dom4j.Element;
import org.dom4j.io.SAXReader;
import proteomics.Types.Dash;
import proteomics.Types.Peptide;
import proteomics.Types.PosMassMap;
import proteomics.Types.VarPtm;

import java.io.InputStream;
import java.text.DecimalFormat;
import java.util.*;

import static proteomics.PIPI.lszDebugScanNum;

public class InferPTM {
    private static int MaxPtmNumInPart;
    public static final byte N_PART = 0;
    public static final byte C_PART = 1;
    public static final byte ANYWHERE = 4;
    public static final byte PEPC = 3;
    public static final byte PEPN = 2;
    public static final byte PROTC = 1;
    public static final byte PROTN = 0;
    public final static DecimalFormat df3 = new DecimalFormat("0.000");
    private final MassTool massTool;
    private final Map<Character, Double> massTable;
    private final Map<Character, Double> fixModMap;
    private final Set<VarPtm> varPtmSet = new HashSet<>();
    private final double minPtmMass;
    private final double maxPtmMass;
    private final double ms2Tol;
    private final Map<Character, List<VarPtm>> finalPtmMap = new HashMap<>(22);
    private Map<Character, Map<Byte, List<VarPtm>>> aaAllVarPtmMap = new HashMap<>(22);
    private final Set<Character> aaCharSet = new HashSet<>(Arrays.asList('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'Y'));

    private final Set<Character> aaWithFixModSet = new HashSet<>();
    private final int g_thread_num;
    private final double timeLimit;
    private final boolean nTermSpecific;
    private final boolean cTermSpecific;

    private final Set<Double> ambiguousMass = new HashSet<>(Arrays.asList(79.957, 114.032, -89.029920));

    public InferPTM(MassTool massTool, Map<Character, Double> fixModMap, Map<String, String> parameterMap) throws Exception {
        nTermSpecific = Integer.valueOf(parameterMap.get("n_specific")) == 1;
        cTermSpecific = Integer.valueOf(parameterMap.get("c_specific")) == 1;
        this.massTool = massTool;
        MaxPtmNumInPart = Integer.valueOf(parameterMap.get("num_sec_ptm"));
        massTable = massTool.getMassTable();
        this.fixModMap = fixModMap;
        for (Character c : fixModMap.keySet()) {
            if (Math.abs(fixModMap.get(c)) > 0.02) {
                aaWithFixModSet.add(c);
            }
        }
        this.minPtmMass = Double.valueOf(parameterMap.get("min_ptm_mass"));
        this.g_thread_num = Integer.valueOf(parameterMap.get("GUROBI_thread_num"));
        this.maxPtmMass = Double.valueOf(parameterMap.get("max_ptm_mass"));
        this.ms2Tol = Double.valueOf(parameterMap.get("ms2_tolerance"));
        this.timeLimit = Double.valueOf(parameterMap.get("GUROBI_time_limit"));

        char[] aaArray = new char[]{'G', 'A', 'S', 'P', 'V', 'T', 'C', 'L', 'N', 'D', 'Q', 'K', 'E', 'M', 'H', 'F', 'R', 'Y', 'W'};
        int n_varPtm = 0;
        for (String k : parameterMap.keySet()) {
            if (!k.startsWith("mod")) continue;
            String v = parameterMap.get(k);
            if (v.startsWith("0.0")) break;
            n_varPtm++;
        }

        Set<Double> enrichedMass = new HashSet<>();
        for (String k : parameterMap.keySet()) {
            if (!k.startsWith("mod")) continue;

            String v = parameterMap.get(k);
            if (v.startsWith("0.0")) break;

            String[] modStr = v.split(",");
            double modMass = Double.valueOf(modStr[0]);
            char modSite = modStr[1].charAt(0);
            byte modPosition = Byte.valueOf(modStr[2]);
            int priority = 1;

            enrichedMass.add(modMass);
            if (modSite == 'M' && modStr[3].contentEquals("Oxidation") && n_varPtm != 0) {
                priority = 0;
            }
            if (modPosition == ANYWHERE) {
                if (Math.abs(fixModMap.get(modSite)) < 0.1) {
                    varPtmSet.add(new VarPtm(modMass, modSite, modPosition, modStr[3], "ByUser", priority));
                }
            } else {
                if (modSite == 'X') {
                    for (char oriAa : aaArray) {
                        if (Math.abs(fixModMap.get(modSite)) < 0.1) {
                            varPtmSet.add(new VarPtm(modMass, oriAa, modPosition, modStr[3], "ByUser", priority));
                        }
                    }
                } else {
                    if (Math.abs(fixModMap.get(modSite)) < 0.1) {
                        varPtmSet.add(new VarPtm(modMass, modSite, modPosition, modStr[3], "ByUser", priority));
                    }
                }
            }
        }

        readModFromUnimod(enrichedMass);

        for (VarPtm varPtm : varPtmSet) {
            if (finalPtmMap.containsKey(varPtm.site)) {
                finalPtmMap.get(varPtm.site).add(varPtm);
            } else {
                List<VarPtm> tempList = new LinkedList<>();
                tempList.add(varPtm);
                finalPtmMap.put(varPtm.site, tempList);
            }
        }

        aaAllVarPtmMap = new HashMap<>(22);
        for (Character aa : finalPtmMap.keySet()) {
            Map<Byte, List<VarPtm>> positionVarPtmMap = new HashMap<>(5);
            for (VarPtm varPtm : finalPtmMap.get(aa)) {
                byte position = varPtm.position;
                if (positionVarPtmMap.containsKey(position)) {
                    positionVarPtmMap.get(position).add(varPtm);
                } else {
                    List<VarPtm> varPtmList = new ArrayList<>(50);
                    varPtmList.add(varPtm);
                    positionVarPtmMap.put(position, varPtmList);
                }
            }
            aaAllVarPtmMap.put(aa, positionVarPtmMap);
        }
        int a = 1;
    }


    public void findBestPtmMIPExtCPWL(int scanNum, GRBEnv env, double totalDeltaMass, int refPos, String partSeq,
                                      double ms1TolAbs, Map<Integer, Integer> absPosYIdMap, Map<Integer, TreeMap<Double, VarPtm>> absPos_MassVarPtm_Map,
                                      int optStartPos,
                                      boolean couldBeProtC, Map<Integer, Integer> yIdMaxAbsPosMap,
                                      TreeMap<Double, Double> plMap, TreeSet<Peptide> cModPepsSet, Map<Integer, List<Byte>> absPos_ptmPositions_Map) {


        Map<Integer, List<Integer>> yIdAllPosesMap = new HashMap<>(absPosYIdMap.values().size());
        for (int pos : absPosYIdMap.keySet()) {
            int yId = absPosYIdMap.get(pos);
            List<Integer> allPoses = yIdAllPosesMap.get(yId);
            if (allPoses != null) {
                allPoses.add(pos);
            } else {
                allPoses = new ArrayList<>();
                allPoses.add(pos);
                yIdAllPosesMap.put(yId, allPoses);
            }
        }
        try {
            GRBModel model = new GRBModel(env);
            GRBLinExpr totalNumsOnPep = new GRBLinExpr();
            GRBLinExpr totalFlexiableMass = new GRBLinExpr();
            Map<Integer, GRBVar> yVarMap = new HashMap<>(yIdAllPosesMap.size());
            for (int yId : absPosYIdMap.values()) {

                GRBVar yVar = model.addVar(0, 1, -0.1, GRB.BINARY, "y" + yId);
                yVarMap.put(yId, yVar);
                double aaMass = 0;
                for (int absPos : yIdAllPosesMap.get(yId)) {
                    aaMass += massTool.getMassTable().get(partSeq.charAt(absPos - refPos));
                }
                totalFlexiableMass.addTerm(aaMass, yVarMap.get(yId));
            }

            Map<Integer, Map<Double, GRBVar>> xVarMap = new HashMap<>();
            Map<Integer, List<Double>> xPosPEPC_MassMap = new HashMap<>();
            int partSeqLen = partSeq.length();
            for (int relPos = 0; relPos < partSeqLen; relPos++) {
                char aa = partSeq.charAt(relPos);
                int absPos = refPos + relPos;

                Map<Double, GRBVar> massXVarMap = new HashMap<>();
                xVarMap.put(absPos, massXVarMap);

                List<Byte> positionsToTry = absPos_ptmPositions_Map.get(absPos);

                GRBLinExpr sumX_leq_1orY = new GRBLinExpr();
                Map<Byte, List<VarPtm>> allVarPtmMap = aaAllVarPtmMap.get(aa);
                Set<Double> usedMassForThisAa = new HashSet<>();
                for (Byte position : positionsToTry) {
                    List<VarPtm> varPtmList = allVarPtmMap.get(position);
                    if (varPtmList == null) continue;
                    double ptmMass;
                    for (VarPtm varPtm : varPtmList) {
                        ptmMass = varPtm.mass;
                        if (massTable.get(partSeq.charAt(relPos)) + ptmMass < ms2Tol) continue;
                        if (usedMassForThisAa.contains(ptmMass)) continue;

                        GRBVar xVar = model.addVar(0, 1, -0.1, GRB.BINARY, "x" + absPos + "_" + ptmMass);
                        massXVarMap.put(ptmMass, xVar);
                        sumX_leq_1orY.addTerm(1, xVar);
                        totalFlexiableMass.addTerm(ptmMass, xVar);
                        totalNumsOnPep.addTerm(1, xVar);
                        usedMassForThisAa.add(ptmMass);
                        if (position == PEPC) {
                            List<Double> PEPC_MassList = xPosPEPC_MassMap.get(absPos);
                            if (PEPC_MassList == null) {
                                PEPC_MassList = new ArrayList<>();
                                PEPC_MassList.add(ptmMass);
                                xPosPEPC_MassMap.put(absPos, PEPC_MassList);
                            } else {
                                PEPC_MassList.add(ptmMass);
                            }
                        }
                    }
                }

                if (sumX_leq_1orY.size() != 0) {
                    if (absPosYIdMap.containsKey(absPos)) {
                        int yId = absPosYIdMap.get(absPos);
                        sumX_leq_1orY.addTerm(-1, yVarMap.get(yId));
                        model.addConstr(sumX_leq_1orY, GRB.LESS_EQUAL, 0, "sumX_leq_Y" + absPos);
                    } else {
                        model.addConstr(sumX_leq_1orY, GRB.LESS_EQUAL, 1, "sumX_leq_1" + absPos);
                    }
                }
            }

            List<Integer> yIdList = new ArrayList<>(yVarMap.keySet());
            yIdList.sort(Comparator.naturalOrder());
            int yIdNum = yIdList.size();
            if (yIdNum >= 2) {
                for (int i = 0; i < yIdNum - 1; i++) {
                    int thisYId = yIdList.get(i);
                    int nextYId = yIdList.get(i + 1);
                    GRBLinExpr yiyj = new GRBLinExpr();
                    yiyj.addTerm(1, yVarMap.get(thisYId));
                    yiyj.addTerm(-1, yVarMap.get(nextYId));
                    model.addConstr(yiyj, GRB.GREATER_EQUAL, 0, "yiyj_" + i);

                    int thisMaxPos = yIdMaxAbsPosMap.get(thisYId);
                    if (!hasFixMod(partSeq.charAt(thisMaxPos - refPos))) {
                        for (double ptmMass : xPosPEPC_MassMap.get(thisMaxPos)) {
                            GRBLinExpr yiyjxi = new GRBLinExpr();
                            yiyjxi.addTerm(1, yVarMap.get(thisYId));
                            yiyjxi.addTerm(-1, yVarMap.get(nextYId));
                            yiyjxi.addTerm(-1, xVarMap.get(thisMaxPos).get(ptmMass));
                            model.addConstr(yiyjxi, GRB.GREATER_EQUAL, 0, "yiyjxi" + i + "_" + ptmMass);
                        }
                    }
                }
                int onlyOneYId = yIdList.get(0);
                if (!hasFixMod(partSeq.charAt(optStartPos - 1 - refPos))) {
                    for (double ptmMass : xPosPEPC_MassMap.get(optStartPos - 1)) {
                        GRBLinExpr yixi_1 = new GRBLinExpr();
                        yixi_1.addTerm(1, yVarMap.get(onlyOneYId));
                        yixi_1.addTerm(1, xVarMap.get(optStartPos - 1).get(ptmMass));
                        model.addConstr(yixi_1, GRB.LESS_EQUAL, 1, "yixi_1" + (optStartPos - 1) + "_" + ptmMass);
                    }
                }
            } else if (yIdNum == 1) {
                int onlyOneYId = yIdList.get(0);
                if (!hasFixMod(partSeq.charAt(optStartPos - 1 - refPos))) {
                    for (double ptmMass : xPosPEPC_MassMap.get(optStartPos - 1)) {
                        GRBLinExpr yixi_1 = new GRBLinExpr();
                        yixi_1.addTerm(1, yVarMap.get(onlyOneYId));
                        yixi_1.addTerm(1, xVarMap.get(optStartPos - 1).get(ptmMass));
                        model.addConstr(yixi_1, GRB.LESS_EQUAL, 1, "yixi_1" + (optStartPos - 1) + "_" + ptmMass);
                    }
                }
            }

            model.addConstr(totalFlexiableMass, GRB.LESS_EQUAL, totalDeltaMass + ms1TolAbs, "constrTotalMassLE");
            model.addConstr(totalFlexiableMass, GRB.GREATER_EQUAL, totalDeltaMass - ms1TolAbs, "constrTotalMassGE");
            model.addConstr(totalNumsOnPep, GRB.LESS_EQUAL, Math.min(partSeq.length(), MaxPtmNumInPart), "totalNumsOnPepConstr");
            int numOfTps = Math.min(partSeqLen - 1, 12);
            for (int tpId = 1; tpId <= numOfTps; tpId++) {
                int absPos = refPos + tpId;
                double minYmz;
                double maxYmz;
                int absPosOfNextStart;
                if (absPosYIdMap.containsKey(absPos)) {
                    int thisYId = absPosYIdMap.get(absPos);
                    absPosOfNextStart = yIdMaxAbsPosMap.get(thisYId) + 1;
                } else {
                    absPosOfNextStart = optStartPos;
                }

                if (absPosYIdMap.containsKey(absPos)) {
                    minYmz = 0;
                } else {
                    int aaNumC = Math.min(MaxPtmNumInPart - 1, absPosOfNextStart - absPos);
                    List<Double> bigPtmMassesOnN = new ArrayList<>(absPos - refPos);
                    for (int tmpPos = refPos; tmpPos < absPos; tmpPos++) {
                        if (absPos_MassVarPtm_Map.containsKey(tmpPos) && absPos_MassVarPtm_Map.get(tmpPos).lastKey() > 0) {
                            bigPtmMassesOnN.add(absPos_MassVarPtm_Map.get(tmpPos).lastKey());
                        } else {
                            bigPtmMassesOnN.add(0d);
                        }
                    }
                    bigPtmMassesOnN.sort(Comparator.reverseOrder());
                    double maxPtmMassOnN = 0;
                    for (int tmpPos = 0; tmpPos < Math.min(MaxPtmNumInPart - aaNumC, bigPtmMassesOnN.size()); ++tmpPos) {
                        maxPtmMassOnN += bigPtmMassesOnN.get(tmpPos);
                    }

                    List<Double> smallPtmMassesHere = new ArrayList<>(absPosOfNextStart - absPos);
                    for (int tmpPos = absPos; tmpPos < absPosOfNextStart; tmpPos++) {
                        if (absPos_MassVarPtm_Map.containsKey(tmpPos) && absPos_MassVarPtm_Map.get(tmpPos).firstKey() < 0) {
                            smallPtmMassesHere.add(absPos_MassVarPtm_Map.get(tmpPos).firstKey());
                        } else {
                            smallPtmMassesHere.add(0d);
                        }
                    }
                    smallPtmMassesHere.sort(Comparator.naturalOrder());
                    double minPtmMassHere = 0;
                    for (int i = 0; i < aaNumC; i++) {
                        minPtmMassHere += smallPtmMassesHere.get(i);
                    }
                    double seqMass_TpId_nextStart = massTool.calResidueMass(partSeq.substring(tpId, absPosOfNextStart - refPos));
                    minYmz = Math.max(0, MassTool.PROTON + massTool.H2O + seqMass_TpId_nextStart + Math.max(totalDeltaMass - maxPtmMassOnN, minPtmMassHere));
                }

                int aaNumC = Math.min(MaxPtmNumInPart - 1, partSeqLen - tpId);
                List<Double> smallPtmMassesC = new ArrayList<>();
                double minPtmMassOnN = 0;
                for (int tmpPos = refPos; tmpPos < absPos; tmpPos++) {
                    if (absPos_MassVarPtm_Map.containsKey(tmpPos) && absPos_MassVarPtm_Map.get(tmpPos).firstKey() < 0) {
                        smallPtmMassesC.add(absPos_MassVarPtm_Map.get(tmpPos).firstKey());
                    } else {
                        smallPtmMassesC.add(0d);
                    }
                }
                smallPtmMassesC.sort(Comparator.naturalOrder());
                for (int tmpPos = 0; tmpPos < Math.min(MaxPtmNumInPart - aaNumC, smallPtmMassesC.size()); tmpPos++) {
                    minPtmMassOnN += smallPtmMassesC.get(tmpPos);
                }

                List<Double> bigPtmMassesHere = new ArrayList<>(partSeqLen - tpId);
                for (int tmpPos = absPos; tmpPos < partSeqLen + refPos; tmpPos++) {
                    if (absPos_MassVarPtm_Map.containsKey(tmpPos) && absPos_MassVarPtm_Map.get(tmpPos).lastKey() > 0) {
                        bigPtmMassesHere.add(absPos_MassVarPtm_Map.get(tmpPos).lastKey());
                    } else {
                        bigPtmMassesHere.add(0d);
                    }
                }
                bigPtmMassesHere.sort(Comparator.reverseOrder());
                double maxPtmMassHere = 0;
                for (int i = 0; i < aaNumC; i++) {
                    maxPtmMassHere += bigPtmMassesHere.get(i);
                }

                maxYmz = Math.min(3000, massTool.H2O + MassTool.PROTON + massTool.calResidueMass(partSeq.substring(tpId)) + Math.min(totalDeltaMass - massTool.calResidueMass(partSeq.substring(absPosOfNextStart - refPos)) - minPtmMassOnN, maxPtmMassHere));
                if (minYmz >= maxYmz) {
                    continue;
                }

                GRBVar yMassVar = model.addVar(minYmz, maxYmz, 0, GRB.CONTINUOUS, "yMassVar_" + tpId);
                GRBLinExpr yMassDef = new GRBLinExpr();
                yMassDef.addTerm(1, yMassVar);
                for (int tmpAbsPos = absPos; tmpAbsPos < refPos + partSeqLen; tmpAbsPos++) {
                    for (double tmpMass : xVarMap.get(tmpAbsPos).keySet()) {
                        yMassDef.addTerm(-tmpMass, xVarMap.get(tmpAbsPos).get(tmpMass));
                    }
                    if (absPosYIdMap.containsKey(tmpAbsPos)) {
                        yMassDef.addTerm(-massTable.get(partSeq.charAt(tmpAbsPos - refPos)), yVarMap.get(absPosYIdMap.get(tmpAbsPos)));
                    } else {
                        yMassDef.addConstant(-massTable.get(partSeq.charAt(tmpAbsPos - refPos)));
                    }
                }
                model.addConstr(yMassDef, GRB.EQUAL, MassTool.PROTON + massTool.H2O, "yMassDef_" + tpId);
                NavigableMap<Double, Double> feasiblePeaks = plMap.subMap(minYmz, true, maxYmz, true);
                int nPeaks = feasiblePeaks.size();
                if (nPeaks != 0) {
                    List<Double> pwX_yIonList = new ArrayList<>(4 * nPeaks + 2);
                    List<Double> pwY_yIonList = new ArrayList<>(4 * nPeaks + 2);
                    generatePWLpts1(feasiblePeaks, pwX_yIonList, pwY_yIonList, minYmz, maxYmz);
                    double[] pwX_yIon = pwX_yIonList.stream().mapToDouble(Double::doubleValue).toArray();
                    double[] pwY_yIon = pwY_yIonList.stream().mapToDouble(Double::doubleValue).toArray();
                    model.setPWLObj(yMassVar, pwX_yIon, pwY_yIon);
                }
            }

            model.set(GRB.IntAttr.ModelSense, GRB.MAXIMIZE);
            model.set(GRB.DoubleParam.TimeLimit, timeLimit);
            model.set(GRB.IntParam.ConcurrentMIP, 4);
            model.set(GRB.IntParam.Threads, g_thread_num);
            model.set(GRB.IntParam.PoolSearchMode, 2);
            model.set(GRB.IntParam.PoolSolutions, 5);
            model.set(GRB.DoubleParam.MIPGap, 1e-1);
            model.set(GRB.IntParam.PreSparsify, 1);

            model.optimize();

            switch (model.get(GRB.IntAttr.Status)) {
                case GRB.OPTIMAL:
                    break;
                case GRB.TIME_LIMIT:
                    break;
                default:
                    model.dispose();
                    return;
            }
            int solCount = model.get(GRB.IntAttr.SolCount);
            if (lszDebugScanNum.contains(scanNum)) {
                int a = 1;
            }
            double optObj = model.get(GRB.DoubleAttr.ObjVal);
            double objValThres;
            if (optObj > 0) {
                objValThres = 0.7 * optObj;
            } else {
                objValThres = optObj - 0.2;
            }
            for (int solNum = 0; solNum < solCount; solNum++) {
                model.set(GRB.IntParam.SolutionNumber, solNum);
                double poolObjVal = model.get(GRB.DoubleAttr.PoolObjVal);
                if (poolObjVal < objValThres) {
                    break;
                }
                Map<Integer, Double> thisSol = new HashMap<>();
                for (int absPos = refPos; absPos < partSeq.length() + refPos; absPos++) {
                    if (absPosYIdMap.containsKey(absPos)) {
                        int yId = absPosYIdMap.get(absPos);
                        GRBVar yVar = yVarMap.get(yId);
                        if ((int) Math.round(yVar.get(GRB.DoubleAttr.Xn)) == 1) {
                            thisSol.put(absPos, -1.0);
                            if (xVarMap.containsKey(absPos)) {
                                Map<Double, GRBVar> ptmIdXVarMap = xVarMap.get(absPos);
                                for (double ptmMass : ptmIdXVarMap.keySet()) {
                                    GRBVar xVar = ptmIdXVarMap.get(ptmMass);
                                    if ((int) Math.round(xVar.get(GRB.DoubleAttr.Xn)) == 1) {
                                        thisSol.put(absPos, ptmMass);
                                        break;
                                    }
                                }
                            }
                        }
                    } else {
                        thisSol.put(absPos, -1.0);
                        if (xVarMap.containsKey(absPos)) {
                            Map<Double, GRBVar> ptmIdXVarMap = xVarMap.get(absPos);
                            for (double ptmMass : ptmIdXVarMap.keySet()) {
                                GRBVar xVar = ptmIdXVarMap.get(ptmMass);
                                if ((int) Math.round(xVar.get(GRB.DoubleAttr.Xn)) == 1) {
                                    thisSol.put(absPos, ptmMass);
                                    break;
                                }
                            }
                        }
                    }
                }
                int trueSeqEndAbsPos = Collections.max(thisSol.keySet());
                Peptide tmpPeptide = new Peptide(partSeq.substring(0, trueSeqEndAbsPos + 1 - refPos), false, massTool);
                PosMassMap posMassMap = new PosMassMap();
                for (int absPos : thisSol.keySet()) {
                    if (thisSol.get(absPos) == -1) continue;
                    posMassMap.put(absPos - refPos, thisSol.get(absPos));
                    tmpPeptide.posVarPtmResMap.put(absPos - refPos, absPos_MassVarPtm_Map.get(absPos).get(thisSol.get(absPos)));
                }
                if (!posMassMap.isEmpty()) tmpPeptide.setVarPTM(posMassMap);
                tmpPeptide.setScore(poolObjVal);
                cModPepsSet.add(tmpPeptide);
            }
            model.dispose();
        } catch (GRBException e) {
            System.out.println("Error code: " + e.getErrorCode() + ". " + e.getMessage());
        }
    }

    public void findBestPtmInGapPWL(int scanNum, GRBEnv env, double totalDeltaMass, int refAbsPos, int rTagPosInProt, String gapSeq,
                                    double ms1TolAbs, Map<Integer, TreeMap<Double, VarPtm>> absPos_MassVarPtm_Map,
                                    TreeMap<Double, Double> plMap, TreeSet<Peptide> midModPepsSet, double bIonRefMass, double yIonRefMass, double precursorMass) {

        try {
            GRBModel model = new GRBModel(env);
            GRBLinExpr totalNumsInGap = new GRBLinExpr();
            GRBLinExpr deltaMass = new GRBLinExpr();

            Map<Integer, Map<Double, GRBVar>> xVarMap = new HashMap<>();
            int gapLen = gapSeq.length();
            for (int relPos = 0; relPos < gapLen; relPos++) {
                char aa = gapSeq.charAt(relPos);
                int absPos = refAbsPos + relPos;

                Map<Double, GRBVar> massXVarMap = new HashMap<>();
                xVarMap.put(absPos, massXVarMap);
                GRBLinExpr sumX_leq_1 = new GRBLinExpr();
                TreeMap<Double, VarPtm> massVarPtmMap = absPos_MassVarPtm_Map.get(absPos);
                if (massVarPtmMap == null) continue;
                for (double ptmMass : massVarPtmMap.keySet()) {
                    if (massTable.get(gapSeq.charAt(relPos)) + ptmMass < ms2Tol) continue;
                    GRBVar xVar = model.addVar(0, 1, -0.1, GRB.BINARY, "x_" + absPos + "_" + ptmMass);
                    massXVarMap.put(ptmMass, xVar);
                    sumX_leq_1.addTerm(1, xVar);
                    deltaMass.addTerm(ptmMass, xVar);
                    totalNumsInGap.addTerm(1, xVar);
                }
                if (sumX_leq_1.size() != 0) {
                    model.addConstr(sumX_leq_1, GRB.LESS_EQUAL, 1, "sumX_leq_1" + absPos);
                }
            }
            model.addConstr(deltaMass, GRB.LESS_EQUAL, totalDeltaMass + ms1TolAbs, "constrTotalMassLE");
            model.addConstr(deltaMass, GRB.GREATER_EQUAL, totalDeltaMass - ms1TolAbs, "constrTotalMassGE");
            model.addConstr(totalNumsInGap, GRB.LESS_EQUAL, Math.min(gapSeq.length(), MaxPtmNumInPart), "totalNumsOnPepConstr");
            int numOfTps = Math.min(gapLen - 1, 12);
            NavigableMap<Double, Double> feasiblePeaks;
            Pair<Double, Double> massRange;
            for (int tpId = 1; tpId <= numOfTps; tpId++) {
                GRBLinExpr bySumMass = new GRBLinExpr();

                int absPos = refAbsPos + tpId;
                massRange = getPotentialYmzRange(rTagPosInProt, absPos, absPos_MassVarPtm_Map, gapSeq, refAbsPos, totalDeltaMass);
                double minYmz = massRange.getFirst() + yIonRefMass;
                double maxYmz = massRange.getSecond() + yIonRefMass;
                if (minYmz >= maxYmz) {
                    continue;
                }
                GRBVar yMassVar = model.addVar(minYmz, maxYmz, 0, GRB.CONTINUOUS, "yMassVar_" + tpId);
                GRBLinExpr yMassDef = new GRBLinExpr();
                yMassDef.addTerm(1, yMassVar);
                for (int tmpAbsPos = absPos; tmpAbsPos < rTagPosInProt; tmpAbsPos++) {
                    for (double tmpMass : xVarMap.get(tmpAbsPos).keySet()) {
                        yMassDef.addTerm(-tmpMass, xVarMap.get(tmpAbsPos).get(tmpMass));
                    }
                    yMassDef.addConstant(-massTable.get(gapSeq.charAt(tmpAbsPos - refAbsPos)));
                }
                model.addConstr(yMassDef, GRB.EQUAL, yIonRefMass, "yMassDef_" + tpId);
                feasiblePeaks = plMap.subMap(minYmz, true, maxYmz, true);
                int nPeaks = feasiblePeaks.size();
                if (nPeaks != 0) {
                    List<Double> pwX_yIonList = new ArrayList<>(4 * nPeaks + 2);
                    List<Double> pwY_yIonList = new ArrayList<>(4 * nPeaks + 2);
                    generatePWLpts1(feasiblePeaks, pwX_yIonList, pwY_yIonList, minYmz, maxYmz);
                    double[] pwX_yIon = pwX_yIonList.stream().mapToDouble(Double::doubleValue).toArray();
                    double[] pwY_yIon = pwY_yIonList.stream().mapToDouble(Double::doubleValue).toArray();
                    model.setPWLObj(yMassVar, pwX_yIon, pwY_yIon);
                }
                bySumMass.addTerm(1, yMassVar);
                massRange = getPotentialBmzRange(rTagPosInProt, absPos, absPos_MassVarPtm_Map, gapSeq, refAbsPos, totalDeltaMass);
                double minBmz = massRange.getFirst() + bIonRefMass;
                double maxBmz = massRange.getSecond() + bIonRefMass;
                if (minBmz >= maxBmz) {
                    continue;
                }
                GRBVar bMassVar = model.addVar(minBmz, maxBmz, 0, GRB.CONTINUOUS, "bMassVar_" + tpId);
                GRBLinExpr bMassDef = new GRBLinExpr();
                bMassDef.addTerm(1, bMassVar);
                for (int tmpAbsPos = refAbsPos; tmpAbsPos < absPos; tmpAbsPos++) {
                    for (double tmpMass : xVarMap.get(tmpAbsPos).keySet()) {
                        bMassDef.addTerm(-tmpMass, xVarMap.get(tmpAbsPos).get(tmpMass));
                    }
                    bMassDef.addConstant(-massTable.get(gapSeq.charAt(tmpAbsPos - refAbsPos)));
                }
                model.addConstr(bMassDef, GRB.EQUAL, bIonRefMass, "bMassDef_" + tpId);
                feasiblePeaks = plMap.subMap(minBmz, true, maxBmz, true);
                nPeaks = feasiblePeaks.size();
                if (nPeaks != 0) {
                    List<Double> pwX_bIonList = new ArrayList<>(4 * nPeaks + 2);
                    List<Double> pwY_bIonList = new ArrayList<>(4 * nPeaks + 2);
                    generatePWLpts1(feasiblePeaks, pwX_bIonList, pwY_bIonList, minBmz, maxBmz);
                    double[] pwX_bIon = pwX_bIonList.stream().mapToDouble(Double::doubleValue).toArray();
                    double[] pwY_bIon = pwY_bIonList.stream().mapToDouble(Double::doubleValue).toArray();
                    model.setPWLObj(bMassVar, pwX_bIon, pwY_bIon);
                }
                bySumMass.addTerm(1, bMassVar);

                model.addConstr(bySumMass, GRB.LESS_EQUAL, precursorMass + 2 * MassTool.PROTON + ms1TolAbs, "bySumMassLEQ_" + tpId);
                model.addConstr(bySumMass, GRB.GREATER_EQUAL, precursorMass + 2 * MassTool.PROTON - ms1TolAbs, "bySumMassGEQ_" + tpId);
            }


            model.set(GRB.IntAttr.ModelSense, GRB.MAXIMIZE);
            model.set(GRB.DoubleParam.TimeLimit, timeLimit);
            model.set(GRB.IntParam.ConcurrentMIP, 4);
            model.set(GRB.IntParam.Threads, g_thread_num);
            model.set(GRB.IntParam.PoolSearchMode, 2);
            model.set(GRB.IntParam.PoolSolutions, 5);
            model.set(GRB.IntParam.PreSparsify, 1);
            model.optimize();
            switch (model.get(GRB.IntAttr.Status)) {
                case GRB.OPTIMAL:
                    break;
                case GRB.TIME_LIMIT:
                    break;
                default:
                    model.dispose();
                    return;
            }
            int solCount = model.get(GRB.IntAttr.SolCount);
            if (lszDebugScanNum.contains(scanNum)) {
                int a = 1;
            }
            double optObj = model.get(GRB.DoubleAttr.ObjVal);
            double objValThres;
            if (optObj > 0) {
                objValThres = 0.7 * optObj;
            } else {
                objValThres = optObj - 0.2;
            }
            for (int solNum = 0; solNum < solCount; solNum++) {
                model.set(GRB.IntParam.SolutionNumber, solNum);
                double poolObjVal = model.get(GRB.DoubleAttr.PoolObjVal);
                if (poolObjVal < objValThres) {
                    break;
                }
                Map<Integer, Double> thisSol = new HashMap<>();
                for (int absPos = refAbsPos; absPos < gapSeq.length() + refAbsPos; absPos++) {
                    thisSol.put(absPos, -1.0);
                    if (xVarMap.containsKey(absPos)) {
                        Map<Double, GRBVar> ptmIdXVarMap = xVarMap.get(absPos);
                        for (double ptmMass : ptmIdXVarMap.keySet()) {
                            GRBVar xVar = ptmIdXVarMap.get(ptmMass);
                            if ((int) Math.round(xVar.get(GRB.DoubleAttr.Xn)) == 1) {
                                thisSol.put(absPos, ptmMass);
                                break;
                            }
                        }
                    }
                }

                Peptide tmpPeptide = new Peptide(gapSeq, false, massTool);
                PosMassMap posMassMap = new PosMassMap();
                for (int absPos : thisSol.keySet()) {
                    if (thisSol.get(absPos) == -1) continue;
                    posMassMap.put(absPos - refAbsPos, thisSol.get(absPos));
                    tmpPeptide.posVarPtmResMap.put(absPos - refAbsPos, absPos_MassVarPtm_Map.get(absPos).get(thisSol.get(absPos)));
                }
                if (!posMassMap.isEmpty()) {
                    tmpPeptide.setVarPTM(posMassMap);
                }
                tmpPeptide.setScore(poolObjVal);
                midModPepsSet.add(tmpPeptide);
            }
            model.dispose();
        } catch (GRBException e) {
            System.out.println("Error code: " + e.getErrorCode() + ". " + e.getMessage());
        }
    }

    private void merge2Dash(Dash lD, Dash rD) {
        if (lD.h < rD.h) {
            lD.y = rD.x;
        } else if (lD.h > rD.h) {
            rD.x = lD.y;
        } else {
            lD.y = rD.y;
            rD.y = rD.x;
        }
    }

    private void generatePWLpts1(NavigableMap<Double, Double> feasiblePeaks, List<Double> pwX_List, List<Double> pwY_List, double minMz, double maxMz) {
        double lastY = 0;
        List<List<Dash>> dashGroupList = new ArrayList<>();
        List<Dash> tmpList = new ArrayList<>();
        for (double mz : feasiblePeaks.keySet()) {
            double thisX = mz - 1.1 * ms2Tol;
            double thisY = mz + 1.1 * ms2Tol;
            if (thisX > lastY && !tmpList.isEmpty()) {
                List<Dash> okList = new ArrayList<>();
                okList.addAll(tmpList);
                dashGroupList.add(okList);
                tmpList.clear();
            }
            tmpList.add(new Dash(thisX, thisY, feasiblePeaks.get(mz)));
            lastY = thisY;
        }
        dashGroupList.add(tmpList);
        for (List<Dash> dashes : dashGroupList) {
            int groupSize = dashes.size();
            Dash lD, rD;
            if (groupSize >= 3) {
                for (int lId = 0; lId < groupSize - 1; lId++) {
                    lD = dashes.get(lId);
                    if (lD.y <= lD.x) continue;
                    for (int rId = lId + 1; rId < groupSize; rId++) {
                        if (lD.y <= lD.x) break;
                        rD = dashes.get(rId);
                        if (rD.y <= rD.x) continue;
                        if (rD.x > lD.y) break;
                        merge2Dash(lD, rD);
                    }
                }
            } else if (groupSize == 2) {
                lD = dashes.get(0);
                rD = dashes.get(1);
                merge2Dash(lD, rD);
            }
        }

        lastY = Math.max(minMz - 400, 0);
        for (List<Dash> dashes : dashGroupList) {
            pwX_List.add(lastY);
            pwY_List.add(0d);
            pwX_List.add(dashes.get(0).x);
            pwY_List.add(0d);
            lastY = dashes.get(dashes.size() - 1).y;
            for (Dash dash : dashes) {
                if (dash.y <= dash.x) continue;
                pwX_List.add(dash.x);
                pwY_List.add(dash.h);
                pwX_List.add(dash.y);
                pwY_List.add(dash.h);
            }
        }
        pwX_List.add(lastY);
        pwY_List.add(0d);
        pwX_List.add(Math.min(lastY + 400, maxMz + 2 * ms2Tol));
        pwY_List.add(0d);
    }

    private Pair<Double, Double> getPotentialYmzRange(int rTagPosInProt, int absPos, Map<Integer, TreeMap<Double, VarPtm>> absPos_MassVarPtm_Map, String gapSeq, int refAbsPos, double totalDeltaMass) {
        double minYmz;
        double maxYmz;
        int cLen = rTagPosInProt - absPos;
        int nLen = absPos - refAbsPos;
        int ptmNumOnC = Math.min(MaxPtmNumInPart - 1, cLen);
        int ptmNumOnN = Math.min(MaxPtmNumInPart - ptmNumOnC, nLen);

        List<Double> allSmallPtmsOnN = new ArrayList<>(nLen);
        List<Double> allBigPtmsOnN = new ArrayList<>(nLen);
        for (int tmpPos = refAbsPos; tmpPos < absPos; tmpPos++) {
            if (absPos_MassVarPtm_Map.containsKey(tmpPos) && absPos_MassVarPtm_Map.get(tmpPos).lastKey() > 0) {
                allBigPtmsOnN.add(absPos_MassVarPtm_Map.get(tmpPos).lastKey());
            } else {
                allBigPtmsOnN.add(0d);
            }
            if (absPos_MassVarPtm_Map.containsKey(tmpPos) && absPos_MassVarPtm_Map.get(tmpPos).firstKey() < 0) {
                allSmallPtmsOnN.add(absPos_MassVarPtm_Map.get(tmpPos).firstKey());
            } else {
                allSmallPtmsOnN.add(0d);
            }
        }
        allBigPtmsOnN.sort(Comparator.reverseOrder());
        allSmallPtmsOnN.sort(Comparator.naturalOrder());
        double maxPtmMassOnN = 0;
        double minPtmMassOnN = 0;
        for (int tmpPos = 0; tmpPos < ptmNumOnN; ++tmpPos) {
            minPtmMassOnN += allSmallPtmsOnN.get(tmpPos);
            maxPtmMassOnN += allBigPtmsOnN.get(tmpPos);
        }

        List<Double> allBigPtmsOnC = new ArrayList<>(cLen);
        List<Double> allSmallPtmsOnC = new ArrayList<>(cLen);
        double sumAaMassOnC = 0;
        for (int tmpPos = absPos; tmpPos < rTagPosInProt; tmpPos++) {
            if (absPos_MassVarPtm_Map.containsKey(tmpPos) && absPos_MassVarPtm_Map.get(tmpPos).firstKey() < 0) {
                allSmallPtmsOnC.add(absPos_MassVarPtm_Map.get(tmpPos).firstKey());
            } else {
                allSmallPtmsOnC.add(0d);
            }
            if (absPos_MassVarPtm_Map.containsKey(tmpPos) && absPos_MassVarPtm_Map.get(tmpPos).lastKey() > 0) {
                allBigPtmsOnC.add(absPos_MassVarPtm_Map.get(tmpPos).lastKey());
            } else {
                allBigPtmsOnC.add(0d);
            }
            sumAaMassOnC += massTable.get(gapSeq.charAt(tmpPos - refAbsPos));
        }
        allSmallPtmsOnC.sort(Comparator.naturalOrder());
        allBigPtmsOnC.sort(Comparator.reverseOrder());
        double minPtmMassOnC = 0;
        double maxPtmMassOnC = 0;
        for (int i = 0; i < ptmNumOnC; i++) {
            minPtmMassOnC += allSmallPtmsOnC.get(i);
            maxPtmMassOnC += allBigPtmsOnC.get(i);
        }
        minYmz = Math.max(0, sumAaMassOnC + Math.max(totalDeltaMass - maxPtmMassOnN, minPtmMassOnC));
        maxYmz = Math.min(3000, sumAaMassOnC + Math.min(totalDeltaMass - minPtmMassOnN, maxPtmMassOnC));
        return new Pair<>(minYmz, maxYmz);
    }

    private Pair<Double, Double> getPotentialBmzRange(int rTagPosInProt, int absPos, Map<Integer, TreeMap<Double, VarPtm>> absPos_MassVarPtm_Map, String gapSeq, int refAbsPos, double totalDeltaMass) {
        double minBmz;
        double maxBmz;
        int nLen = absPos - refAbsPos;
        int cLen = rTagPosInProt - absPos;
        int ptmNumOnN = Math.min(MaxPtmNumInPart - 1, nLen);
        int ptmNumOnC = Math.min(MaxPtmNumInPart - ptmNumOnN, cLen);
        List<Double> allSmallPtmsOnN = new ArrayList<>(nLen);
        List<Double> allBigPtmsOnN = new ArrayList<>(nLen);
        double sumAaMassOnN = 0;
        for (int tmpPos = refAbsPos; tmpPos < absPos; tmpPos++) {
            if (absPos_MassVarPtm_Map.containsKey(tmpPos) && absPos_MassVarPtm_Map.get(tmpPos).lastKey() > 0) {
                allBigPtmsOnN.add(absPos_MassVarPtm_Map.get(tmpPos).lastKey());
            } else {
                allBigPtmsOnN.add(0d);
            }
            if (absPos_MassVarPtm_Map.containsKey(tmpPos) && absPos_MassVarPtm_Map.get(tmpPos).firstKey() < 0) {
                allSmallPtmsOnN.add(absPos_MassVarPtm_Map.get(tmpPos).firstKey());
            } else {
                allSmallPtmsOnN.add(0d);
            }
            sumAaMassOnN += massTable.get(gapSeq.charAt(tmpPos - refAbsPos));
        }
        allBigPtmsOnN.sort(Comparator.reverseOrder());
        allSmallPtmsOnN.sort(Comparator.naturalOrder());
        double maxPtmMassOnN = 0;
        double minPtmMassOnN = 0;
        for (int tmpPos = 0; tmpPos < ptmNumOnN; ++tmpPos) {
            minPtmMassOnN += allSmallPtmsOnN.get(tmpPos);
            maxPtmMassOnN += allBigPtmsOnN.get(tmpPos);
        }

        List<Double> allBigPtmsOnC = new ArrayList<>(cLen);
        List<Double> allSmallPtmsOnC = new ArrayList<>(cLen);
        for (int tmpPos = absPos; tmpPos < rTagPosInProt; tmpPos++) {
            if (absPos_MassVarPtm_Map.containsKey(tmpPos) && absPos_MassVarPtm_Map.get(tmpPos).firstKey() < 0) {
                allSmallPtmsOnC.add(absPos_MassVarPtm_Map.get(tmpPos).firstKey());
            } else {
                allSmallPtmsOnC.add(0d);
            }
            if (absPos_MassVarPtm_Map.containsKey(tmpPos) && absPos_MassVarPtm_Map.get(tmpPos).lastKey() > 0) {
                allBigPtmsOnC.add(absPos_MassVarPtm_Map.get(tmpPos).lastKey());
            } else {
                allBigPtmsOnC.add(0d);
            }
        }
        allSmallPtmsOnC.sort(Comparator.naturalOrder());
        allBigPtmsOnC.sort(Comparator.reverseOrder());
        double minPtmMassOnC = 0;
        double maxPtmMassOnC = 0;
        for (int i = 0; i < ptmNumOnC; i++) {
            minPtmMassOnC += allSmallPtmsOnC.get(i);
            maxPtmMassOnC += allBigPtmsOnC.get(i);
        }

        minBmz = Math.max(0, sumAaMassOnN + Math.max(totalDeltaMass - maxPtmMassOnN, minPtmMassOnC));
        maxBmz = Math.min(3000, sumAaMassOnN + Math.min(totalDeltaMass - minPtmMassOnN, maxPtmMassOnC));
        return new Pair<>(minBmz, maxBmz);
    }

    public void findBestPtmMIPExtNPWL(int scanNum, GRBEnv env, double totalDeltaMass, int refPos, String partSeq,
                                      double ms1TolAbs, Map<Integer, Integer> absPosYIdMap,
                                      Map<Integer, TreeMap<Double, VarPtm>> absPos_MassVarPtm_Map, int optEndPosP1,
                                      boolean couldBeProtN, Map<Integer, Integer> yIdMinAbsPosMap,
                                      String protSeq, TreeMap<Double, Double> plMap,
                                      TreeSet<Peptide> nModPepsSet,
                                      Map<Integer, List<Byte>> absPos_ptmPositions_Map) {

        Map<Integer, List<Integer>> yIdAllPosesMap = new HashMap<>();
        for (int pos : absPosYIdMap.keySet()) {
            int yId = absPosYIdMap.get(pos);
            List<Integer> allPoses = yIdAllPosesMap.get(yId);
            if (allPoses != null) {
                allPoses.add(pos);
            } else {
                allPoses = new ArrayList<>();
                allPoses.add(pos);
                yIdAllPosesMap.put(yId, allPoses);
            }
        }

        try {
            GRBModel model = new GRBModel(env);
            model.set(GRB.IntAttr.ModelSense, GRB.MAXIMIZE);

            GRBLinExpr totalNumsOnPep = new GRBLinExpr();
            GRBLinExpr totalFlexiableMass = new GRBLinExpr();
            int partSeqLen = partSeq.length();
            Map<Integer, GRBVar> yVarMap = new HashMap<>(yIdAllPosesMap.size());
            for (int yId : yIdAllPosesMap.keySet()) {
                GRBVar yVar = model.addVar(0, 1, -0.1, GRB.BINARY, "y" + yId);
                yVarMap.put(yId, yVar);
                double aaMass = 0;
                for (int absPos : yIdAllPosesMap.get(yId)) {
                    aaMass += massTool.getMassTable().get(partSeq.charAt(absPos + partSeqLen - refPos));
                }
                totalFlexiableMass.addTerm(aaMass, yVarMap.get(yId));
            }

            Map<Integer, Map<Double, GRBVar>> xVarMap = new HashMap<>();
            Map<Integer, List<Double>> xPosPEPN_MassMap = new HashMap<>();
            for (int relPos = 0; relPos < partSeqLen; relPos++) {
                char aa = partSeq.charAt(relPos);
                int absPos = refPos + relPos - partSeqLen;


                GRBLinExpr sumX_leq_1orY = new GRBLinExpr();
                Map<Byte, List<VarPtm>> allVarPtmMap = aaAllVarPtmMap.get(aa);
                Map<Double, GRBVar> massXVarMap = new HashMap<>();
                xVarMap.put(absPos, massXVarMap);
                Set<Double> usedMassForThisAa = new HashSet<>();
                List<Byte> positionsToTry = absPos_ptmPositions_Map.get(absPos);
                for (Byte position : positionsToTry) {
                    List<VarPtm> varPtmList = allVarPtmMap.get(position);
                    if (varPtmList == null) continue;

                    double ptmMass;
                    for (VarPtm varPtm : varPtmList) {
                        ptmMass = varPtm.mass;
                        if (massTable.get(partSeq.charAt(relPos)) + ptmMass < ms2Tol) continue;
                        if (usedMassForThisAa.contains(ptmMass)) continue;

                        GRBVar xVar = model.addVar(0, 1, -0.1, GRB.BINARY, "x" + absPos + "_" + ptmMass);
                        massXVarMap.put(ptmMass, xVar);
                        sumX_leq_1orY.addTerm(1, xVar);
                        totalFlexiableMass.addTerm(ptmMass, xVar);
                        totalNumsOnPep.addTerm(1, xVar);
                        usedMassForThisAa.add(ptmMass);
                        if (position == PEPN) {
                            List<Double> PEPN_MassList = xPosPEPN_MassMap.get(absPos);
                            if (PEPN_MassList == null) {
                                PEPN_MassList = new ArrayList<>();
                                PEPN_MassList.add(ptmMass);
                                xPosPEPN_MassMap.put(absPos, PEPN_MassList);
                            } else {
                                PEPN_MassList.add(ptmMass);
                            }
                        }
                    }
                }
                if (sumX_leq_1orY.size() != 0) {
                    if (absPosYIdMap.containsKey(absPos)) {
                        int yId = absPosYIdMap.get(absPos);
                        sumX_leq_1orY.addTerm(-1, yVarMap.get(yId));
                        model.addConstr(sumX_leq_1orY, GRB.LESS_EQUAL, 0, "sumX_leq_Y" + absPos);
                    } else {
                        model.addConstr(sumX_leq_1orY, GRB.LESS_EQUAL, 1, "sumX_leq_1" + absPos);
                    }
                }
            }

            List<Integer> yIdList = new ArrayList<>(yVarMap.keySet());
            yIdList.sort(Comparator.naturalOrder());
            int yIdNum = yIdList.size();
            if (yIdNum >= 2) {
                for (int i = 0; i < yIdNum - 1; i++) {
                    int thisYId = yIdList.get(i);
                    int nextYId = yIdList.get(i + 1);
                    GRBLinExpr yiyj = new GRBLinExpr();
                    yiyj.addTerm(1, yVarMap.get(thisYId));
                    yiyj.addTerm(-1, yVarMap.get(nextYId));
                    model.addConstr(yiyj, GRB.GREATER_EQUAL, 0, "yiyj_" + i);

                    int thisMinPos = yIdMinAbsPosMap.get(thisYId);
                    if (!hasFixMod(partSeq.charAt(thisMinPos + partSeqLen - refPos))) {
                        for (double ptmMass : xPosPEPN_MassMap.get(thisMinPos)) {
                            GRBLinExpr yiyjxi = new GRBLinExpr();
                            yiyjxi.addTerm(1, yVarMap.get(thisYId));
                            yiyjxi.addTerm(-1, yVarMap.get(nextYId));
                            yiyjxi.addTerm(-1, xVarMap.get(thisMinPos).get(ptmMass));
                            model.addConstr(yiyjxi, GRB.GREATER_EQUAL, 0, "yiyjxi" + i + "_" + ptmMass);
                        }
                    }
                }
                int onlyOneYId = yIdList.get(0);
                if (!hasFixMod(partSeq.charAt(optEndPosP1 + partSeqLen - refPos))) {
                    for (double ptmMass : xPosPEPN_MassMap.get(optEndPosP1)) {
                        GRBLinExpr yixi_1 = new GRBLinExpr();
                        yixi_1.addTerm(1, yVarMap.get(onlyOneYId));
                        yixi_1.addTerm(1, xVarMap.get(optEndPosP1).get(ptmMass));
                        model.addConstr(yixi_1, GRB.LESS_EQUAL, 1, "yixi_1" + (optEndPosP1) + "_" + ptmMass);
                    }
                }
            } else if (yIdNum == 1) {
                int onlyOneYId = yIdList.get(0);
                if (!hasFixMod(partSeq.charAt(optEndPosP1 + partSeqLen - refPos))) {
                    for (double ptmMass : xPosPEPN_MassMap.get(optEndPosP1)) {
                        GRBLinExpr yixi_1 = new GRBLinExpr();
                        yixi_1.addTerm(1, yVarMap.get(onlyOneYId));
                        yixi_1.addTerm(1, xVarMap.get(optEndPosP1).get(ptmMass));
                        model.addConstr(yixi_1, GRB.LESS_EQUAL, 1, "yixi_1" + (optEndPosP1) + "_" + ptmMass);
                    }
                }
            }

            model.addConstr(totalFlexiableMass, GRB.LESS_EQUAL, totalDeltaMass + ms1TolAbs, "constrTotalMassLE");
            model.addConstr(totalFlexiableMass, GRB.GREATER_EQUAL, totalDeltaMass - ms1TolAbs, "constrTotalMassGE");
            model.addConstr(totalNumsOnPep, GRB.LESS_EQUAL, Math.min(partSeqLen, MaxPtmNumInPart), "totalNumsOnPep");
            int numOfTps = Math.min(partSeqLen - 1, 12);
            for (int tpId = 1; tpId <= numOfTps; tpId++) {
                int absPos = refPos + tpId - partSeqLen;
                double minBmz;
                double maxBmz;
                int absPosOfNextEnd;
                if (absPosYIdMap.containsKey(absPos)) {
                    int thisYId = absPosYIdMap.get(absPos);
                    absPosOfNextEnd = yIdMinAbsPosMap.get(thisYId) - 1;
                } else {
                    absPosOfNextEnd = optEndPosP1 - 1;
                }

                if (absPosYIdMap.containsKey(absPos) || absPosYIdMap.containsKey(absPos - 1)) {
                    minBmz = 0;
                } else {
                    int aaNumN = Math.min(MaxPtmNumInPart - 1, absPos - absPosOfNextEnd - 1);
                    List<Double> bigPtmMassesOnC = new ArrayList<>(refPos - absPos - 1);
                    for (int tmpPos = absPos; tmpPos < refPos; tmpPos++) {
                        if (absPos_MassVarPtm_Map.containsKey(tmpPos) && absPos_MassVarPtm_Map.get(tmpPos).lastKey() > 0) {
                            bigPtmMassesOnC.add(absPos_MassVarPtm_Map.get(tmpPos).lastKey());
                        } else {
                            bigPtmMassesOnC.add(0d);
                        }
                    }
                    bigPtmMassesOnC.sort(Comparator.reverseOrder());
                    double maxPtmMassOnC = 0;
                    for (int tmpPos = 0; tmpPos < Math.min(MaxPtmNumInPart - aaNumN, bigPtmMassesOnC.size()); ++tmpPos) {
                        maxPtmMassOnC += bigPtmMassesOnC.get(tmpPos);
                    }

                    List<Double> smallPtmMassesHere = new ArrayList<>(absPos - absPosOfNextEnd);
                    for (int tmpPos = absPosOfNextEnd + 1; tmpPos < absPos; tmpPos++) {
                        if (absPos_MassVarPtm_Map.containsKey(tmpPos) && absPos_MassVarPtm_Map.get(tmpPos).firstKey() < 0) {
                            smallPtmMassesHere.add(absPos_MassVarPtm_Map.get(tmpPos).firstKey());
                        } else {
                            smallPtmMassesHere.add(0d);
                        }
                    }
                    smallPtmMassesHere.sort(Comparator.naturalOrder());
                    double minPtmMassHere = 0;
                    for (int i = 0; i < aaNumN; i++) {
                        minPtmMassHere += smallPtmMassesHere.get(i);
                    }
                    double seqMass_nextEnd_TpId = massTool.calResidueMass(partSeq.substring(absPosOfNextEnd + 1 - refPos + partSeqLen, tpId));
                    minBmz = Math.max(0, MassTool.PROTON + seqMass_nextEnd_TpId + Math.max(totalDeltaMass - maxPtmMassOnC, minPtmMassHere));
                }

                int aaNumN = Math.min(MaxPtmNumInPart - 1, tpId + 1);
                List<Double> smallPtmMassesC = new ArrayList<>();
                double minPtmMassOnC = 0;
                for (int tmpPos = absPos + 1; tmpPos < refPos; tmpPos++) {
                    if (absPos_MassVarPtm_Map.containsKey(tmpPos) && absPos_MassVarPtm_Map.get(tmpPos).firstKey() < 0) {
                        smallPtmMassesC.add(absPos_MassVarPtm_Map.get(tmpPos).firstKey());
                    } else {
                        smallPtmMassesC.add(0d);
                    }
                }
                smallPtmMassesC.sort(Comparator.naturalOrder());
                for (int tmpPos = 0; tmpPos < Math.min(MaxPtmNumInPart - aaNumN, smallPtmMassesC.size()); tmpPos++) {
                    minPtmMassOnC += smallPtmMassesC.get(tmpPos);
                }
                List<Double> bigPtmMassesHere = new ArrayList<>(tpId + 1);
                for (int tmpPos = refPos - partSeqLen; tmpPos < absPos + 1; tmpPos++) {
                    if (absPos_MassVarPtm_Map.containsKey(tmpPos) && absPos_MassVarPtm_Map.get(tmpPos).lastKey() > 0) {
                        bigPtmMassesHere.add(absPos_MassVarPtm_Map.get(tmpPos).lastKey());
                    } else {
                        bigPtmMassesHere.add(0d);
                    }
                }
                bigPtmMassesHere.sort(Comparator.reverseOrder());
                double maxPtmMassHere = 0;
                for (int i = 0; i < aaNumN; i++) {
                    maxPtmMassHere += bigPtmMassesHere.get(i);
                }
                maxBmz = Math.min(3000, MassTool.PROTON + massTool.calResidueMass(partSeq.substring(0, tpId + 1)) + Math.min(totalDeltaMass - massTool.calResidueMass(partSeq.substring(0, absPosOfNextEnd + 1 - refPos + partSeqLen)) - minPtmMassOnC, maxPtmMassHere));
                if (minBmz > maxBmz) {
                    continue;
                }

                GRBVar bMassVar = model.addVar(minBmz, maxBmz, 0, GRB.CONTINUOUS, "bMassVar_" + tpId);
                GRBLinExpr bMassDef = new GRBLinExpr();
                bMassDef.addTerm(1, bMassVar);
                for (int tmpAbsPos = refPos - partSeqLen; tmpAbsPos < absPos; tmpAbsPos++) {
                    for (double tmpMass : xVarMap.get(tmpAbsPos).keySet()) {
                        bMassDef.addTerm(-tmpMass, xVarMap.get(tmpAbsPos).get(tmpMass));
                    }
                    if (absPosYIdMap.containsKey(tmpAbsPos)) {
                        bMassDef.addTerm(-massTable.get(partSeq.charAt(tmpAbsPos - refPos + partSeqLen)), yVarMap.get(absPosYIdMap.get(tmpAbsPos)));
                    } else {
                        bMassDef.addConstant(-massTable.get(partSeq.charAt(tmpAbsPos - refPos + partSeqLen)));
                    }
                }
                model.addConstr(bMassDef, GRB.EQUAL, MassTool.PROTON, "bMassDef_" + tpId);
                NavigableMap<Double, Double> feasiblePeaks = plMap.subMap(minBmz, true, maxBmz, true);
                int nPeaks = feasiblePeaks.size();
                if (nPeaks != 0) {
                    List<Double> pwX_bIonList = new ArrayList<>(4 * nPeaks + 2);
                    List<Double> pwY_bIonList = new ArrayList<>(4 * nPeaks + 2);
                    generatePWLpts1(feasiblePeaks, pwX_bIonList, pwY_bIonList, minBmz, maxBmz);
                    double[] pwX_bIon = pwX_bIonList.stream().mapToDouble(Double::doubleValue).toArray();
                    double[] pwY_bIon = pwY_bIonList.stream().mapToDouble(Double::doubleValue).toArray();
                    model.setPWLObj(bMassVar, pwX_bIon, pwY_bIon);
                }
            }

            model.set(GRB.DoubleParam.TimeLimit, timeLimit);
            model.set(GRB.IntParam.ConcurrentMIP, 4);
            model.set(GRB.IntParam.Threads, g_thread_num);
            model.set(GRB.IntParam.PreSparsify, 1);
            model.set(GRB.IntParam.PoolSearchMode, 2);
            model.set(GRB.IntParam.PoolSolutions, 5);


            model.optimize();
            switch (model.get(GRB.IntAttr.Status)) {
                case GRB.OPTIMAL:
                    break;
                case GRB.TIME_LIMIT:
                    break;
                default:
                    model.dispose();
                    return;
            }
            int solCount = model.get(GRB.IntAttr.SolCount);
            double optObj = model.get(GRB.DoubleAttr.ObjVal);
            double objValThres;
            if (optObj > 0) {
                objValThres = 0.7 * optObj;
            } else {
                objValThres = optObj - 0.2;
            }
            for (int solNum = 0; solNum < solCount; solNum++) {
                model.set(GRB.IntParam.SolutionNumber, solNum);
                double poolObjVal = model.get(GRB.DoubleAttr.PoolObjVal);
                if (model.get(GRB.DoubleAttr.PoolObjVal) < objValThres) {
                    break;
                }
                Map<Integer, Double> thisSol = new HashMap<>();
                for (int absPos = refPos - partSeqLen; absPos < refPos; absPos++) {
                    if (absPosYIdMap.containsKey(absPos)) {
                        int yId = absPosYIdMap.get(absPos);
                        GRBVar yVar = yVarMap.get(yId);
                        if ((int) Math.round(yVar.get(GRB.DoubleAttr.Xn)) == 1) {
                            thisSol.put(absPos, -1.0);
                            if (xVarMap.containsKey(absPos)) {
                                Map<Double, GRBVar> ptmIdXVarMap = xVarMap.get(absPos);
                                for (double ptmMass : ptmIdXVarMap.keySet()) {
                                    GRBVar xVar = ptmIdXVarMap.get(ptmMass);
                                    if ((int) Math.round(xVar.get(GRB.DoubleAttr.Xn)) == 1) {
                                        thisSol.put(absPos, ptmMass);
                                        break;
                                    }
                                }
                            }
                        }
                    } else {
                        thisSol.put(absPos, -1.0);
                        if (xVarMap.containsKey(absPos)) {
                            Map<Double, GRBVar> ptmIdXVarMap = xVarMap.get(absPos);
                            for (double ptmMass : ptmIdXVarMap.keySet()) {
                                GRBVar xVar = ptmIdXVarMap.get(ptmMass);
                                if ((int) Math.round(xVar.get(GRB.DoubleAttr.Xn)) == 1) {
                                    thisSol.put(absPos, ptmMass);
                                    break;
                                }
                            }
                        }
                    }
                }

                int trueSeqStartAbsPos = Collections.min(thisSol.keySet());
                Peptide tmpPeptide = new Peptide(partSeq.substring(trueSeqStartAbsPos - refPos + partSeqLen), false, massTool);
                PosMassMap posMassMap = new PosMassMap();
                for (int absPos : thisSol.keySet()) {
                    if (thisSol.get(absPos) == -1) continue;
                    posMassMap.put(absPos - trueSeqStartAbsPos, thisSol.get(absPos));
                    tmpPeptide.posVarPtmResMap.put(absPos - trueSeqStartAbsPos, absPos_MassVarPtm_Map.get(absPos).get(thisSol.get(absPos)));
                }
                if (!posMassMap.isEmpty()) tmpPeptide.setVarPTM(posMassMap);
                tmpPeptide.setScore(poolObjVal);
                nModPepsSet.add(tmpPeptide);
                if (java.lang.management.ManagementFactory.getRuntimeMXBean().getInputArguments().toString().indexOf("jdwp") >= 0) {
                    double[][] ions = tmpPeptide.getIonMatrixNow();
                    int a = 1;
                }
            }
            model.dispose();
        } catch (GRBException e) {
            System.out.println(scanNum + "," + partSeq + " , Error code: " + e.getErrorCode() + ". " + e.getMessage());
        }
    }

    public void updateIonMatrix(double[][] ionMatrix, double cutMass, byte ncPart) {
        if (ncPart == N_PART) {
            for (int i = 0; i < ionMatrix[1].length; i++) {
                ionMatrix[1][i] += cutMass;
            }
        } else {
            for (int i = 0; i < ionMatrix[0].length; i++) {
                ionMatrix[0][i] += cutMass;
            }
        }
    }

    public double getMinPtmMass() {
        return minPtmMass;
    }

    public double getMaxPtmMass() {
        return maxPtmMass;
    }

    private void readModFromUnimod(Set<Double> enrichedMass) throws Exception {
        SAXReader reader = new SAXReader();
        InputStream inputStream = getClass().getClassLoader().getResourceAsStream("unimod.xml");
        Document document = reader.read(inputStream);
        Element rootElement = document.getRootElement();
        Iterator<Element> rootIter = rootElement.elementIterator();

        while (rootIter.hasNext()) {
            Element rootElem = rootIter.next();
            if (!rootElem.getName().contentEquals("modifications")) continue;

            Iterator<Element> modIter = rootElem.elementIterator();

            while (modIter.hasNext()) {
                Element modElem = modIter.next();

                String name = modElem.attributeValue("title");
                if (name.contentEquals("Diethylphosphothione")) {
                    int a = 1;
                }
                double mass = Double.valueOf(modElem.element("delta").attributeValue("mono_mass"));
                if (mass < minPtmMass || mass > maxPtmMass) continue;
                boolean shouldKeep = true;
                for (double tmpmass : ambiguousMass) {
                    if (Math.abs(mass - tmpmass) < 0.005) {
                        shouldKeep = false;
                    }
                }
                if (!shouldKeep) continue;

                shouldKeep = true;
                for (double tmpmass : enrichedMass) {
                    if (Math.abs(mass - tmpmass) < ms2Tol) {
                        shouldKeep = false;
                    }
                }
                if (!shouldKeep) continue;

                for (Element spec : modElem.elements("specificity")) {
                    String classification = spec.attributeValue("classification");
                    if (classification.contains("glycos") || classification.contains("Other")) {
                        continue;
                    }

                    String siteStr = spec.attributeValue("site");
                    String positionStr = spec.attributeValue("position");
                    if (classification.contentEquals("Isotopic label")
                            && !(name.contentEquals("Propionyl") && siteStr.contentEquals("K"))
                            && !(name.contentEquals("Succinyl") && siteStr.contentEquals("K"))
                            && !name.contentEquals("Dimethyl:2H(4)13C(2)")) {
                        continue;
                    }
                    byte position = 0;
                    switch (positionStr) {
                        case "Protein N-term":
                            position = PROTN;
                            break;
                        case "Protein C-term":
                            position = PROTC;
                            break;
                        case "Any N-term":
                            position = PEPN;
                            break;
                        case "Any C-term":
                            position = PEPC;
                            break;
                        case "Anywhere":
                            position = ANYWHERE;
                            break;
                    }
                    if (siteStr.contentEquals("N-term") || siteStr.contentEquals("C-term")) {
                        for (char site : aaCharSet) {
                            if (aaWithFixModSet.contains(site) && siteStr.contentEquals("C-term")) {
                                continue;
                            }

                            VarPtm temp = new VarPtm(mass, site, position, name, classification, 0);
                            if (finalPtmMap.containsKey(site)) {
                                finalPtmMap.get(site).add(temp);
                            } else {
                                List<VarPtm> varPtmSet = new LinkedList<>();
                                varPtmSet.add(temp);
                                finalPtmMap.put(site, varPtmSet);
                            }
                        }
                    } else {
                        char site = siteStr.charAt(0);
                        if (aaWithFixModSet.contains(site) && (position == 1 || position == 3 || position == 4)) {
                            continue;
                        }
                        VarPtm temp = new VarPtm(mass, site, position, name, classification, 0);
                        if (finalPtmMap.containsKey(site)) {
                            finalPtmMap.get(site).add(temp);
                        } else {
                            List<VarPtm> varPtmSet = new LinkedList<>();
                            varPtmSet.add(temp);
                            finalPtmMap.put(site, varPtmSet);
                        }
                    }
                }
            }
        }
    }


    private boolean hasFixMod(char aa) {
        return fixModMap.containsKey(aa) && Math.abs(fixModMap.get(aa)) > 0.1;
    }


    public void prepareInfoCTerm(int scanNum, String partSeq,
                                 Map<Integer, TreeMap<Double, VarPtm>> pos_MassVarPtm_Map,
                                 final Map<Integer, Integer> yIdMaxAbsPosMap,
                                 final boolean couldBeProtC,
                                 final int optStartPos,
                                 final int startRefPos,
                                 Map<Integer, List<Byte>> absPos_ptmPositions_Map) {

        int partSeqLen = partSeq.length();
        char aa;
        int absPos;
        for (int relPos = 0; relPos < partSeqLen; ++relPos) {
            aa = partSeq.charAt(relPos);
            absPos = relPos + startRefPos;
            List<Byte> positionsToTry = new ArrayList<>(3);
            absPos_ptmPositions_Map.put(absPos, positionsToTry);
            if (aaWithFixModSet.contains(aa)) continue;

            positionsToTry.add(ANYWHERE);
            if (yIdMaxAbsPosMap.containsValue(absPos)) {
                positionsToTry.add(PEPC);
            }
            if (absPos == optStartPos - 1) {
                if ((aa == 'K' || aa == 'R') || !cTermSpecific) {
                    positionsToTry.add(PEPC);
                }
            }
            if (couldBeProtC && relPos == partSeqLen - 1) {
                positionsToTry.add(PROTC);
            }

            if (aaAllVarPtmMap.containsKey(aa)) {
                Map<Byte, List<VarPtm>> allVarPtmMap = aaAllVarPtmMap.get(aa);
                Map<String, VarPtm> dstMap = new HashMap<>(128);
                for (Byte position : positionsToTry) {
                    if (!allVarPtmMap.containsKey(position)) continue;
                    for (VarPtm varPtm : allVarPtmMap.get(position)) {
                        double mass = varPtm.mass;
                        if (massTable.get(partSeq.charAt(relPos)) + mass < ms2Tol) continue;
                        String varPtmStr = varPtm.getStr();

                        VarPtm oldVarPtm = dstMap.get(varPtmStr);
                        if (oldVarPtm != null) {
                            if (varPtm.priority > oldVarPtm.priority) {
                                dstMap.put(varPtmStr, varPtm);
                            }
                        } else {
                            dstMap.put(varPtmStr, varPtm);
                        }
                    }
                }
                if (!dstMap.isEmpty()) {
                    TreeMap<Double, VarPtm> massVarPtmMap = new TreeMap<>();
                    for (VarPtm varPtm : dstMap.values()) {
                        massVarPtmMap.put(varPtm.mass, varPtm);
                    }
                    pos_MassVarPtm_Map.put(absPos, massVarPtmMap);
                }
            }
        }
    }

    public void prepareInfoMid(int scanNum, String partSeq, Map<Integer, TreeMap<Double, VarPtm>> pos_MassVarPtm_Map, final int startRefPos) {

        int partSeqLen = partSeq.length();
        char aa;
        int absPos;
        for (int relPos = 0; relPos < partSeqLen; ++relPos) {
            aa = partSeq.charAt(relPos);
            absPos = relPos + startRefPos;
            if (aaWithFixModSet.contains(aa)) continue;
            if (aaAllVarPtmMap.containsKey(aa)) {
                Map<Byte, List<VarPtm>> allVarPtmMap = aaAllVarPtmMap.get(aa);
                Map<String, VarPtm> dstMap = new HashMap<>(128);
                if (!allVarPtmMap.containsKey(ANYWHERE)) continue;
                for (VarPtm varPtm : allVarPtmMap.get(ANYWHERE)) {
                    double mass = varPtm.mass;
                    if (massTable.get(partSeq.charAt(relPos)) + mass < ms2Tol) continue;
                    String varPtmStr = varPtm.getStr();

                    VarPtm oldVarPtm = dstMap.get(varPtmStr);
                    if (oldVarPtm != null) {
                        if (varPtm.priority > oldVarPtm.priority) {
                            dstMap.put(varPtmStr, varPtm);
                        }
                    } else {
                        dstMap.put(varPtmStr, varPtm);
                    }
                }
                if (!dstMap.isEmpty()) {
                    TreeMap<Double, VarPtm> massVarPtmMap = new TreeMap<>();
                    for (VarPtm varPtm : dstMap.values()) {
                        massVarPtmMap.put(varPtm.mass, varPtm);
                    }
                    pos_MassVarPtm_Map.put(absPos, massVarPtmMap);
                }
            }
        }
    }

    public void findPossible1Ptm(int scanNum, String partSeq, Map<Integer, TreeMap<Double, VarPtm>> pos_MassVarPtm_Map, int startRefPos, TreeSet<Peptide> modPepsSet, double deltaMass, double ms1TolAbs) {

        int partSeqLen = partSeq.length();
        int absPos;
        for (int relPos = 0; relPos < partSeqLen; ++relPos) {
            absPos = relPos + startRefPos;
            TreeMap<Double, VarPtm> massVarPtmMap = pos_MassVarPtm_Map.get(absPos);
            if (massVarPtmMap == null) continue;
            for (double mass : massVarPtmMap.subMap(deltaMass - ms1TolAbs - ms2Tol, deltaMass + ms1TolAbs + ms2Tol).keySet()) {
                if (massTable.get(partSeq.charAt(relPos)) + mass < ms2Tol) continue;
                VarPtm varPtm = massVarPtmMap.get(mass);

                if ((varPtm.position == 2 || varPtm.position == 0) && relPos != 0) continue;
                if ((varPtm.position == 3 || varPtm.position == 1) && relPos != partSeqLen - 1) continue;
                Peptide tmpPeptide = new Peptide(partSeq, false, massTool);
                PosMassMap posMassMap = new PosMassMap();
                posMassMap.put(relPos, mass);
                tmpPeptide.posVarPtmResMap.put(relPos, varPtm);
                tmpPeptide.setVarPTM(posMassMap);
                tmpPeptide.setScore(-0.1 * posMassMap.size());
                modPepsSet.add(tmpPeptide);
            }
        }
    }

    public void prepareInfoNTerm(int scanNum, String partSeq,
                                 Map<Integer, TreeMap<Double, VarPtm>> absPos_MassVarPtm_Map,
                                 Map<Integer, Integer> yIdMinAbsPosMap,
                                 boolean couldBeProtN,
                                 int optEndPosP1,
                                 int endRefPos,
                                 String protSeq,
                                 Map<Integer, List<Byte>> absPos_ptmPositions_Map) {

        int partSeqLen = partSeq.length();
        char aa;
        int absPos;
        for (int relPos = 0; relPos < partSeqLen; ++relPos) {
            aa = partSeq.charAt(relPos);
            absPos = relPos + endRefPos - partSeqLen;
            List<Byte> positionsToTry = new ArrayList<>(3);
            absPos_ptmPositions_Map.put(absPos, positionsToTry);
            if (aaWithFixModSet.contains(aa) && absPos != 0) continue;
            positionsToTry.add(ANYWHERE);
            if (yIdMinAbsPosMap.containsValue(absPos)) {
                positionsToTry.add(PEPN);
            }
            if (absPos == optEndPosP1) {
                positionsToTry.add(PEPN);
            }
            if (couldBeProtN && relPos == 0) {
                positionsToTry.add(PROTN);
            }

            if (aaAllVarPtmMap.containsKey(aa)) {
                Map<Byte, List<VarPtm>> allVarPtmMap = aaAllVarPtmMap.get(aa);
                Map<String, VarPtm> dstMap = new HashMap<>();
                for (Byte position : positionsToTry) {
                    if (!allVarPtmMap.containsKey(position)) continue;
                    for (VarPtm varPtm : allVarPtmMap.get(position)) {
                        double mass = varPtm.mass;
                        if (massTable.get(partSeq.charAt(relPos)) + mass < ms2Tol) continue;
                        String varPtmStr = varPtm.getStr();

                        VarPtm oldVarPtm = dstMap.get(varPtmStr);
                        if (oldVarPtm != null) {
                            if (varPtm.priority > oldVarPtm.priority) {
                                dstMap.put(varPtmStr, varPtm);
                            }
                        } else {
                            dstMap.put(varPtmStr, varPtm);
                        }

                    }
                }

                if (!dstMap.isEmpty()) {
                    TreeMap<Double, VarPtm> massVarPtmMap = new TreeMap<>();
                    for (VarPtm varPtm : dstMap.values()) {
                        massVarPtmMap.put(varPtm.mass, varPtm);
                    }
                    absPos_MassVarPtm_Map.put(absPos, massVarPtmMap);
                }
            }
        }
    }

}
