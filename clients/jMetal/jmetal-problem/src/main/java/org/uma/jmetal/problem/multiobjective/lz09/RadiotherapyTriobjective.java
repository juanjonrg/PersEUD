package org.uma.jmetal.problem.multiobjective.lz09;

import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.io.IOException;

import org.uma.jmetal.problem.doubleproblem.impl.AbstractDoubleProblem;
import org.uma.jmetal.solution.Solution;
import org.uma.jmetal.solution.doublesolution.DoubleSolution;
import org.uma.jmetal.util.JMetalException;
import org.zeromq.SocketType;
import org.zeromq.ZContext;
import org.zeromq.ZMQ;
import org.zeromq.ZMQ.Poller;
import org.zeromq.ZMQ.Socket;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;


import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

/**
 * Class representing problem LZ09F9
 */
/** Class representing problem Radiotherapy */

@SuppressWarnings("serial")
public class RadiotherapyTriobjective extends AbstractDoubleProblem {

  private final static int REQUEST_TIMEOUT = 500000; // msecs, (> 1000!)
  private final static int REQUEST_RETRIES = 3; // Before we abandon
  private final static String SERVER_ENDPOINT = "tcp://0.0.0.0:3557";

  private Socket client;
  private final ZContext context;


  /** Constructor */
  public RadiotherapyTriobjective() {
    setNumberOfVariables(15);
    setNumberOfObjectives(2);
    setNumberOfConstraints(0);
    setName("RadiotherapyTriobjective");

    List<Double> lowerLimit = new ArrayList<>(getNumberOfVariables());
    List<Double> upperLimit = new ArrayList<>(getNumberOfVariables());

    // Limits Slinianka P and L /////////////////////////
    for (int i = 0; i < 2; i++) {
      // dosis prescrita EUD0
      lowerLimit.add(0.5);
      upperLimit.add(26.0);
      // a
      lowerLimit.add(1.0);
      upperLimit.add(100.0);

      // pentalty
      lowerLimit.add(1.0);
      upperLimit.add(100.0);
    }

    // Spina /////////////////////////
    for (int i = 0; i < 1; i++) {
      // dosis prescrita EUD0
      lowerLimit.add(0.5);
      upperLimit.add(50.0);
      // a
      lowerLimit.add(1.0);
      upperLimit.add(50.0);
      // pentalty
      lowerLimit.add(1.0);
      upperLimit.add(100.0);
    }

    // TUMORS /////////////////////////////////////
    for (int i = 0; i < 3; i++) {
      // a
      lowerLimit.add(-100.0);
      upperLimit.add(-1.0);
      // penalty
      lowerLimit.add(1.0);
      upperLimit.add(100.0);
    }

    setVariableBounds(lowerLimit, upperLimit);

    // Socket to talk to server
    context = new ZContext();
    System.out.println("Connecting to server");

    client = context.createSocket(SocketType.REQ);
    client.connect(SERVER_ENDPOINT);

  }
  
  /** Evaluate() method */
  public void evaluate(List<DoubleSolution> solutionList) {

    Map<String, Object> dataset = new HashMap<>();
    dataset.put("Router", SERVER_ENDPOINT);
    dataset.put("Patient", "5");
    dataset.put("Study", "TriObj-SGSC-20231027");
    List<Map<String, Object>> plans = new ArrayList<>();

    // Iterate each solution/plan
    for(int i = 0; i < solutionList.size(); i++) {
      DoubleSolution solution = solutionList.get(i);
     
      Map<String, Object> plan = new HashMap<>();
      plan.put("ID", i);

      Map<String, Map<String, Double>> variables = new HashMap<>();
      variables.put("Salivary Gland L", createParameterMap(solution.getVariable(0), solution.getVariable(1), solution.getVariable(2)));
      variables.put("Salivary Gland R", createParameterMap(solution.getVariable(3), solution.getVariable(4), solution.getVariable(5)));
      variables.put("Spinal Cord +3mm", createParameterMap(solution.getVariable(6), solution.getVariable(7), solution.getVariable(8)));
      variables.put("PTV 54", createParameterMap(solution.getVariable(9), solution.getVariable(10)));
      variables.put("PTV 60", createParameterMap(solution.getVariable(11), solution.getVariable(12)));
      variables.put("PTV 66", createParameterMap(solution.getVariable(13), solution.getVariable(14)));
      
      plan.put("Variables", variables);

      Map<String, Double> objectives_old = new HashMap<>();
      objectives_old.put("f0", solution.getObjective(0));
      objectives_old.put("f1", solution.getObjective(1));
      plan.put("Objectives", objectives_old);

      plans.add(plan);
    }

    dataset.put("Plans", plans);

    // Convert the dataset to a JSON string
    Gson gson = new GsonBuilder().create();
    String json = gson.toJson(dataset);
    //json = json.replaceAll("\r", "").replaceAll("\n", "").replaceAll("\t", "").replaceAll(" ", "");
    System.out.println(json);
    // Send json with all the information of the solution
    client.send(json.getBytes(ZMQ.CHARSET), 0);

    // Get the reply
    byte[] reply = client.recv(0);
    String string_reply = new String(reply, ZMQ.CHARSET);
    System.out.println(string_reply);
    List<String> objectives_new = new ArrayList<String>();
    try {
      objectives_new = readAndSplitObjectives(string_reply);
      
    } catch (IOException e) {
      e.printStackTrace();
    }
    
    for (int i = 0; i < solutionList.size(); i++) {
      String[] parts = objectives_new.get(i).split(",");
      try {
        solutionList.get(i).setObjective(0, Double.parseDouble(parts[1]));
        solutionList.get(i).setObjective(1, Double.parseDouble(parts[3]));

      } catch (NumberFormatException e) {
        // Handle the exception as per your requirement
        System.err.println("Error parsing double: " + e.getMessage());
        solutionList.get(i).setObjective(0, 10000.0);
        solutionList.get(i).setObjective(1, 10000.0);
     }
   
    }
  
  }

  public static List<String> readAndSplitObjectives(String jsonData) throws IOException {
        List<String> objectivesList = new ArrayList<>();

        ObjectMapper objectMapper = new ObjectMapper();
        JsonNode root = objectMapper.readTree(jsonData);

        if (root.isArray()) {
            for (JsonNode element : root) {
                JsonNode objective1Node = element.get("f0");
                JsonNode objective2Node = element.get("f1");

                if (objective1Node != null && objective2Node != null) {
                    String objective1 = objective1Node.asText();
                    String objective2 = objective2Node.asText();

                    String combinedObjective = "F 0," + objective1 + ", F 1: ," + objective2;
                    objectivesList.add(combinedObjective);
                }
            }
        }

        return objectivesList;
    }

  private static Map<String, Double> createParameterMap(double a, double n) {
    Map<String, Double> map = new HashMap<>();
    map.put("a", a);
    map.put("n", n);
    return map;
  }

  private static Map<String, Double> createParameterMap(double eud0, double a, double n) {
    Map<String, Double> map = createParameterMap(a, n);
    map.put("EUD_0", eud0);
    return map;

  }

  @Override
  public void evaluate(DoubleSolution solution) {
    // TODO Auto-generated method stub

  }
}
