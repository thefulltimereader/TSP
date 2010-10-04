import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.Scanner;

public class Tsp {
  private Node[] nodes;
  private Double[][] costMatrix;
  private Double[][] pheromone; 
  private int total;
  private double learningRate;
  private double randNot;
  private boolean goBack;
  public Tsp(int length, boolean back){
    total = length;
    nodes = new Node[total];
    costMatrix = new Double[total][total];
    pheromone = new Double[total][total];
    learningRate = 0.5;
    randNot = ((double)new Random().nextInt(100))/100;
    goBack = back;
  }
  public void solve(){
    System.out.println("Start..");
    buildMatrix();
    initPheromone();
    antColony(4, 3, 120);
    
  }
  private void antColony(int samples, int beamWidth, double numChildren){
    ArrayList<Integer> bestSofar = new ArrayList<Integer>();
    ArrayList<Integer> restartBest = new ArrayList<Integer>();
    ArrayList<Integer> iterationBest = new ArrayList<Integer>();
    double convgFactor = 0;
    boolean bs_update = false;//goes true when alg reaches convergence
    long t = System.currentTimeMillis();
    long end = (long) (t + 1.3*60*1000); //20 sec..2 min: 120000.
    int i = 0;
    while(System.currentTimeMillis() < end){
      iterationBest = beamSearch(beamWidth, numChildren, samples);
      //iterationBest = localSearch(iterationBest);
      if(totalCost(iterationBest)< totalCost(restartBest)){
        restartBest = iterationBest;
      }
      if(totalCost(iterationBest)< totalCost(bestSofar)){
        bestSofar = iterationBest;
      }
      
      convgFactor = computeConvergenceFactor();
      
      if(bs_update && convgFactor>0.99){
        initPheromone();
        restartBest = new ArrayList<Integer>();
        bs_update = false;
      }
      else{
        if(convgFactor>0.99) bs_update = true;
        applyPheromoneUpdate(convgFactor, bs_update, iterationBest, restartBest, bestSofar);
      } 
      System.out.println("iteration: "+i);
      i++;
    }
    if(goBack){
      bestSofar.add(0);
    }
    System.out.println("Final solution: "+ bestSofar);
    System.out.println("with cost " + totalCost(bestSofar));
    showResults(bestSofar);
    
    
  }
  /**
   * the three computed solution (tour) determines how to update the pheromone
   * the influence of each depends on the convg factor
   * each pheromone value T_ij is updated as:
   * t_ij = t_ij + p*(gam_ij - t_ij), 
   * where gam_ij = iterBK*iterationBest_ij + restBestK*restartBest_ij + bestFarK * bestFar_ij
   * according to table 1:
   * bs_update == false:
   * if cf <0.4 iterBK = 1, restBestK = 0, bestFarK = 0
   * if cf <0.6 iterBK = 2/3, restBestK = 1/3, bestFarK = 0
   * if cf<0.8 iterBk = 1/3, restBestK = 2/3, bestFarK = 0
   * if cf<1 iterBk = 0, restBestK = 1, bestFarK =0
   * bs_update == true:

   * @param convgFactor
   * @param bsUpdate
   * @param iterationBest
   * @param restartBest
   * @param bestSofar
   */
  private void applyPheromoneUpdate(double convgFactor, boolean bsUpdate,
      ArrayList<Integer> iterationBest, ArrayList<Integer> restartBest,
      ArrayList<Integer> bestSofar) {
    Double[][] gamma = getGamma(iterationBest, restartBest, bestSofar, bsUpdate, convgFactor);
    for(int i =0; i<total; i++){
      for(int j=0; j<total; j++){
        pheromone[i][j] = pheromone[i][j] + learningRate*(gamma[i][j]-pheromone[i][j]);
        //System.out.println("new pheromone at i:"+ i+ " j:" + j+ " is " + pheromone[i][j]);
        if(pheromone[i][j]>0.999) pheromone[i][j] = 0.999;
        if(pheromone[i][j]<0.001) pheromone[i][j] = 0.001;
      }
    }
    
  }
  private Double[][] getGamma(ArrayList<Integer> iterationBest,
      ArrayList<Integer> restartBest, ArrayList<Integer> bestSofar, 
      boolean bsUpdate, double convgFactor) {
    Double[][] gamma = new Double[total][total];
    if(bsUpdate){
      for(int i=0; i<total; i++){
        for(int j= 0; j<total; j++){
          gamma[i][j] = getBestBool(bestSofar, i, j)*1;
        }
      }
   }
    else{
      for(int i=0; i<total; i++){
        for(int j= 0; j<total; j++){
          if(convgFactor<0.4){
            gamma[i][j] = getBestBool(iterationBest, i, j)*1;
          }
          else if(convgFactor<0.6){
            gamma[i][j] = getBestBool(iterationBest, i, j)*(2/3) +
            getBestBool(restartBest, i, j)*(1/3);
          }
          else if(convgFactor<0.8){
            gamma[i][j] = getBestBool(iterationBest, i, j)*(1/3) +
            getBestBool(restartBest, i, j)*(2/3);
          }
          else{
            gamma[i][j] = getBestBool(restartBest, i, j)*1;
          }
        }
      }
    }
    return gamma;
  }
  /** iterBk = restBestK = 0 = bestFarK = 1
  * 
  * xxBest_ij is 1 if j is visited after i in solution xxBest
  */
  private double getBestBool(ArrayList<Integer> bestOf, int i, int j){
    if( bestOf.indexOf(i) <  bestOf.indexOf(j)) return 1; 
    else return 0;
  }
  /**
   * Computes the convg factor which is a function of the current pheromone values
   * @param pheromone
   * @return
   */
  private double computeConvergenceFactor() {
    double sum = 0;
    double maxP = max(pheromone);
    double minP = min(pheromone);
    double sumOfAll = 0;
    for(int i=0; i<total;i++){
      for(int j=0; j<total;j++){
        sum+=maxOfTwo(maxP-pheromone[i][j], pheromone[i][j]-minP);
        //System.out.println(maxP+ " minP: " +minP + " ij" + pheromone[i][j]);
        sumOfAll += pheromone[i][j];
      }
    }
    if(maxP-minP==0) sum = 0.5;
    else sum = sum / (sumOfAll*(maxP-minP));
    double result = 2*(sum-0.5);
    if(result>1){
      result = result/5;
    }
    return result;
  }
  private ArrayList<Integer> localSearch(ArrayList<Integer> path){
    
    for(int k = 1; k<total-1; k++){
      Random ranGen = new Random();
      int rand = ranGen.nextInt(total-3) + 1; 
      ArrayList<Integer> newPath = swap(path, k, rand);
      if(totalCost(newPath) < totalCost(path)){
        path = newPath;
        System.out.println("Swap nodes k:" + path.get(k) + " and "+ path.get(rand));
      }
    }
    return path;
  }
  @SuppressWarnings("unchecked")
  private ArrayList<Integer> swap(ArrayList<Integer> path, int k, int r) {
    ArrayList<Integer> exp = (ArrayList<Integer>) path.clone();
    int temp = exp.get(k);
    int index = exp.indexOf(temp);
    exp.remove(index);
    exp.add(index, exp.get(r));
    int indexR =exp.indexOf(exp.get(r)); 
    exp.remove(indexR);
    exp.add(indexR, temp);   
    return exp;
  }
 
  private ArrayList<Integer> beamSearch(int beamWidth, double numChildren, int samples){
    //this carries each path (partial sols)
    ArrayList<ArrayList<Integer>>partialSols = new ArrayList<ArrayList<Integer>>();
    ArrayList<Integer> start = new ArrayList<Integer>(); start.add(0);
    partialSols.add(start);
    ArrayList<CostAndPath> bt1 = new ArrayList<CostAndPath>();
    
    for(int t=0; t<total-1; t++){
      //C:=C(B_t)
      ArrayList<ArrayList<Integer>> poss = constPossibleExtension(partialSols);
      int bound = (int) Math.min(poss.size(), Math.floor(numChildren*beamWidth));
      for(int k = 1; k<bound; k++){
        //<P,j> = ChooseFrom(C)
        ArrayList<Integer> best = chooseFrom(poss, t);
        //System.out.println("best is: " + best);
        poss.remove(best);
        CostAndPath candp = new CostAndPath(totalCost(best), best);
        bt1.add(candp);
      }
      partialSols.clear();
      partialSols = reduce(bt1, beamWidth);
      bt1.clear();
    }
   // System.out.println(partialSols.get(0));
    return partialSols.get(0);
  }
  
  private ArrayList<ArrayList<Integer>> reduce(ArrayList<CostAndPath> befor, int width){
    Collections.sort(befor);
    //cut it down
    if(width+1<befor.size()) befor.subList(width+1, befor.size()).clear();
    //System.out.println(befor);
    ArrayList<ArrayList<Integer>> nextBatch = new ArrayList<ArrayList<Integer>>();
    for(CostAndPath cap: befor){
      nextBatch.add(cap.getPath());
    }
    return nextBatch;
  }

  private ArrayList<ArrayList<Integer>> constPossibleExtension(ArrayList<ArrayList<Integer>> paths){
    ArrayList<ArrayList<Integer>> newPaths = new ArrayList<ArrayList<Integer>>();
    for(int i=0; i<total; i++){
      for(int j = 0; j<paths.size(); j++){
        ArrayList<Integer> aPoss = new ArrayList<Integer>();
        aPoss.addAll(paths.get(j));
        if(!aPoss.contains(i)){
          aPoss.add(i);
          newPaths.add(aPoss);
        }
      }
    }
    return newPaths;
  }
  
  private ArrayList<Integer> chooseFrom(ArrayList<ArrayList<Integer>> poss, int index){
    double max = 0;
    int best=0;
    ArrayList<Integer> bestPath = null;
    
    Random randGen = new Random();
    int rand = randGen.nextInt(100); 
    double p = ((double)rand)/100;
    if(p<randNot){
      //Deteministic
      for(Iterator<ArrayList<Integer>> itr = poss.iterator(); itr.hasNext();){
        ArrayList<Integer> aPath = itr.next();
        Integer j = aPath.get(aPath.size()-1);//get the last city visited
        //System.out.println(aPath);
        Integer i = aPath.get(index);
        if(j!=i){
          //System.out.println("cost of "+ i+" and "+j+" is " + costMatrix[i][j]);
          double pheromoneCost = pheromone[i][j]/costMatrix[i][j];
          if (pheromoneCost > max){
            max = pheromoneCost;
            best = j;
            bestPath = aPath;
          }
        }
      }
      return bestPath;
    }
    //do stochastic
    else{
      double[] probs = buildProbArray(poss, index);
      return getPathOnProbability(poss, probs, p);

    }

  }
  private ArrayList<Integer> getPathOnProbability(ArrayList<ArrayList<Integer>> poss,
      double[] probs, double rand){
    
    for(int k =0; k<probs.length; k++){
      rand -=probs[k];
     // System.out.println("rand: "+rand);  
      if(rand<=0) {
        return poss.get(k);
      }
    }
    
    return null;
  }
  private double[] buildProbArray(ArrayList<ArrayList<Integer>> poss, int index){
    double[] probs = new double[poss.size()];
    
    double sum = 0;
    for(ArrayList<Integer> aPath: poss){
      int j = aPath.get(aPath.size()-1);//get the last city visited
      int i = aPath.get(index);// index is which the city before the last city
      sum += pheromoneCost(i, j);
    }
    double t = 0;
    for(int k=0; k<poss.size(); k++){
      ArrayList<Integer> aPath = poss.get(k);
      int j = aPath.get(aPath.size()-1);//get the last city visited
      int i = aPath.get(index);// index is which the city before the last city
      double numerator= pheromoneCost(i,j);
      probs[k] = numerator/sum;
      t +=probs[k];
    }
   
    return probs;
    
  }
  private double pheromoneCost(int i, int j){
    return pheromone[i][j]/costMatrix[i][j];
  }
  public void read(String name, int len) throws FileNotFoundException{ 
    Scanner scanner = null;
    try {
      System.out.println("Reading "+name);
      scanner = new Scanner(new File(name));
      while (scanner.hasNextLine()) {
        int index = scanner.nextInt();
        int x = scanner.nextInt();
        int y = scanner.nextInt();
        int z = scanner.nextInt();
        Node n = new Node(x,y,z);
        //System.out.println(n+" at " + (index-1));
        nodes[index-1] = n;
      }
    } finally {
      scanner.close();
    }
  }
  

  private double getDistance(Node from, Node to) {
    if(from!=to){
      double sum = Math.pow(to.getX() - from.getX(), 2)
      + Math.pow(to.getY() - from.getY(), 2)
      + Math.pow(to.getZ() - from.getZ(), 2);
      return Math.sqrt(sum);
    }
    return Double.MAX_VALUE;
    
  }
  private double totalCost(List<Integer> path){
    if(path.isEmpty()){
      return Double.MAX_VALUE;
    }
    double cost = 0;
    int start = path.get(0);
    //System.out.println(path);
    for(int i=1; i<path.size(); i++){
      if(start!=path.get(i)) cost+=costMatrix[start][path.get(i)];
      start = path.get(i);
    }
    return cost;
  }
  public void showResults(List<Integer> path){
    int before = 0;
    double sum = 0;
    System.out.println(path.get(0)+" "+nodes[path.get(0)].getX() + " " + nodes[path.get(0)].getY()+
        " "+ nodes[path.get(0)].getZ() + " dist: "+costMatrix[before][path.get(0)]);
    for(int i=1; i<path.size(); i++){
      if(before!=path.get(i)){
      System.out.println(path.get(i)+" "+nodes[path.get(i)].getX() + " " + nodes[path.get(i)].getY()+
          " "+ nodes[path.get(i)].getZ() + " dist: "+costMatrix[before][path.get(i)]);
      sum = sum+ costMatrix[before][path.get(i)];
      }
      before = path.get(i);
      
    }
    System.out.println("Cost is: " + sum);
    checkValid(path);
  }
  private void checkValid(List<Integer> path){
    if(path.size()!=total){ System.err.println("bad!");} ;
    for(int i = 0; i<total; i++){
      if(!path.contains(i)){ System.err.println("incomprehensive "+ i +  " missing");}
    }
  }
  private void buildMatrix(){
    for(int i=0; i<total; i++){
      for(int j=0;j<total;j++){
        costMatrix[i][j] = getDistance(nodes[i], nodes[j]);
      }
    }
  }
  private void initPheromone(){
    for(int i=0; i<total; i++){
      for(int j=0;j<total;j++){
        pheromone[i][j] = 0.5;
      }
    }
  }
  private double max(Double[][] mat){
    double max = -1;
    for(int i=1; i<mat.length;i++){
      for(int j= 1; j<mat.length;j++){
        if(mat[i][j] > max) max = mat[i][j];
      }
    }
    return max;
  }
  private double min(Double[][] mat){
    double min = Double.MAX_VALUE;
    for(int i=1; i<mat.length;i++){
      for(int j= 1; j<mat.length;j++){
        if(mat[i][j] < min) min = mat[i][j];
      }
    }
    return min;
  }
  private double maxOfTwo(double a, double b){
    if (a > b) return a;
    return b;
  }
  /** Private Node class **/
  class Node {
    private int x;
    private int y;
    private int z;

    public Node(int x, int y, int z) {
      this.x = x;
      this.y = y;
      this.z = z;
    }

    public int getX() {
      return x;
    }

    public int getY() {
      return y;
    }

    public int getZ() {
      return z;
    }
    
    public String toString(){
      return "x:"+x+" y:"+y+" z:"+z;
    }

    @Override
    public boolean equals(Object obj) {
      return ((Node)obj).getX()==this.x && ((Node)obj).getY()==this.y &&
      ((Node)obj).getY()==this.y;
    }
  }
  class CostAndPath implements Comparable<CostAndPath>{
    private Double cost;
    private ArrayList<Integer> path;
    public CostAndPath(Double cost, ArrayList<Integer> path){
      this.cost = cost;
      this.path = path;
    }
    
    public Double getCost() {
      return cost;
    }

    public ArrayList<Integer> getPath() {
      return path;
    }

    public int compareTo(CostAndPath comp){
      return this.cost.compareTo(comp.cost);
    }
    public String toString(){
      return "cost: " + cost + " path: " + path.toString();
    }
  }
  /**
   * Main
   * @param args
   */
  public static void main(String[] args) {
    String fileName = args[0];
    int len = Integer.valueOf(args[1]);
    String flag = args[2];
    boolean back = false;
    if(flag.equals("yes")) back = true;
    Tsp tspSolver = new Tsp(len, back);
    try {
      tspSolver.read(fileName, len);
     
      tspSolver.solve();
      
    } catch (FileNotFoundException e) {
      System.err.println("File DNE, stop.");
      e.printStackTrace();
    }
  }
}
