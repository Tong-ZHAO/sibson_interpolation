# Sibson Interpolation

Common techniques of interpolation sometimes lead to regularation problems. Sibson interpolation helps us avoid this danger by using a refined interpolation.

## Requirement

* g++ (your version should support c++ 11)
* cgal
* python 2.7

## Function

* Sibson interpolation with points chosen randomly
* Given a interpolation error, choose points one by one randomly
* Given a interpolation error, choose points hierarchically
* Given the number of points, minimise the interpolation error

## Run

```
sh run.sh
cd build
./sibson n eps mode input.bmp output.bmp
```

## At last

For more details, please read the report in the folder **report**
