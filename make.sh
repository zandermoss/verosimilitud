if  [ -a ../Verosim.so ] ; then rm ; rm ../Verosim.so; fi
if  [ -a ../libverosimilitud.so ] ; then rm ; ../libverosimilitud.so ; fi
if  [ -a pyclient/Verosim.so ] ; then rm ; rm pyclient/Verosim.so ; fi
if  [ -a pyclient/Verosim.cpp ] ; then rm ; rm pyclient/Verosim.cpp ; fi
if  [ -a libverosimilitud.so ] ; then rm; rm libverosimilitud.so ; fi
if  [ -d libverosimilitud.so.dSYM ] ; then rm -r; rm -r libverosimilitud.so.dSym ; fi
make
cp libverosimilitud.so ..
cd pyclient
sh pycompile.sh
cp Verosim.so ../..
cd ../..
sudo install_name_tool -change libverosimilitud.so "/Users/marjon/Dropbox (MIT)/work/IceCube/decay/neutrino_decay/libverosimilitud.so" "/Users/marjon/Dropbox (MIT)/work/IceCube/decay/neutrino_decay/Verosim.so"
