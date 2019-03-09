About
-----
Source code of the paper [A Region Search Evolutionary Algorithm for Many-Objective Optimization]  
Code also can be seen from (https://www.researchgate.net/profile/Yongqi_Liu/publications)

If you use the code or find it helpful, please cite the following paper:
```
@inproceedings{Yongqi2019,
    title={A Region Search Evolutionary Algorithm for Many-Objective Optimization},
    author={Yongqi Liu, Hui Qin, Zhendong Zhang, Liqiang Yao, Chao Wang, Li Mo, Shuo Ouyang ,Jie Li},
    booktitle={Information Sciences},
    year={2019}
}
```

## The code is build based on jMetal
**jMetal** is an object-oriented Java-based framework for multi-objective optimization with metaheuristics.
The Web page of the project is: [http://jmetal.github.io/jMetal/](http://jmetal.github.io/jMetal/). Former jMetal versions can be found in [SourceForge](http://jmetal.sourceforge.net). 

## Main code of RSEA
RSEA: /jmetal-algorithm/src/main/java/org/uma/jmetal/algorithm/multiobjective/moead/RSEA.java;  
MOEA/D-RS: /jmetal-algorithm/src/main/java/org/uma/jmetal/algorithm/multiobjective/moead/MOEAD_RS.java.

## How to run
For test: /jmetal-mine/src/main/java/org/mine/exec/RSEARunner.java;  
To run whole result (all test problems, each problem run 20 times): /jmetal-mine/src/main/java/org/mine/exec/RSEARunnerRunResult.java.

