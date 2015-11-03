#!/usr/bin/env python
import sys

__author__ = 'hayssam'
import simplejson


#expects two json file as input, load them, and merge their keys


with open(sys.argv[1],"r") as f:
	json_object1=simplejson.load(f)
with open(sys.argv[2],"r") as f:
	json_object2=simplejson.load(f)
json_object1.update(json_object2)
print simplejson.dumps(json_object1)




