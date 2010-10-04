import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Scanner;

public class Tsp {
  private Node[] nodes;
  private Double[][] costMatrix;
  private Double[] pheromone; 
  private int total;
  
  public Tsp(int length){
    total = length;
    nodes = new Node[total];
    costMatrix = new Double[total][total];
    pheromone = new Double[total];
  }
  public void solve(){
    buildMatrix();
    initPheromone();
    beamSearch(5, 8.5, 4);
    antColony(4, 5, 8.5, 0.3);
    
  }
  private void antColony(int samples, int beamWidth, double numChildren, 
      double rand){
    Integer[] bestSoFar = new Integer[100];
    Integer[] restartBest = new Integer[100];
    double convgFactor = 0;
    boolean bs_update = false;//goes true when alg reaches convergence
    
  }
  private void beamSearch(int beamWidth, double numChildren, int samples){
    //this carries each path (partial sols)
    ArrayList<ArrayList<Integer>>partialSols = new ArrayList<ArrayList<Integer>>();
    ArrayList<Integer> start = new ArrayList<Integer>(); start.add(0);
    partialSols.add(start);
    ArrayList<CostAndPath> bt1 = new ArrayList<CostAndPath>();
    
    for(int t=0; t<total; t++){
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
      //System.out.println(partialSols);
    }
  }
  
  private ArrayList<ArrayList<Integer>> reduce(ArrayList<CostAndPath> befor, int width){
    Collections.sort(befor);
    //cut it down
    if(width+1<befor.size()) befor.subList(width+1, befor.size()).clear();
    System.out.println(befor);
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
    for(Iterator<ArrayList<Integer>> itr = poss.iterator(); itr.hasNext();){
      ArrayList<Integer> aPath = itr.next();
      //System.out.println(aPath);
      Integer j = aPath.get(aPath.size()-1);//get the last city visited
      Integer i = aPath.get(index);
      if(j!=i){
        //System.out.println("cost of "+ i+" and "+j+" is " + costMatrix[i][j]);
        double pheromoneCost = 0.5/costMatrix[i][j];
        if (pheromoneCost > max){
          max = pheromoneCost;
          best = j;
          bestPath = aPath;
        }
      }
    }
    return bestPath;
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
    return -1.0;
    
  }
  private double totalCost(List<Integer> path){
    double cost = 0;
    int start = path.get(0);
    for(int i=1; i<path.size(); i++){
      cost+=costMatrix[start][path.get(i)];
      start = i;
    }
    return cost;
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
      pheromone[i] = 0.5;
    }
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
    Tsp tspSolver = new Tsp(len);
    try {
      tspSolver.read(fileName, len);
     
      tspSolver.solve();
      
    } catch (FileNotFoundException e) {
      System.err.println("File DNE, stop.");
      e.printStackTrace();
    }
  }
}
