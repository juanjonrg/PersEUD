# PersEUD (jMetal + ZMQ + Distributed Gradient Descent)
Requirements: 
```
maven/3.6.3
openjdk/17.0.5_8
```
Compile with:

```
cd clients/jMetal
mvn package
```
And then run with:
```
export CLASSPATH=[PATHTO]/clients/jMetal/jmetal-core/target/jmetal-core-5.10-jar-with-dependencies.jar:[PATHTO]/clients/jMetal/jmetal-problem/target/jmetal-problem-5.10-jar-with-dependencies.jar:[PATHTO]/clients/jMetal/jmetal-example/target/jmetal-example-5.10-jar-with-dependencies.jar:[PATHTO]/clients/jMetal/jmetal-problem/target/jmetal-problem-5.10-jar-with-dependencies.jar
java org.uma.jmetal.example.multiobjective.moead.MOEADRunner
```

## Running the router

```
cd routers
./router.py
```

## Running the multicore worker

Compile with:

```
cd workers/multicore
. env.sh
make
```
And then run with:
```
./run.sh
```

