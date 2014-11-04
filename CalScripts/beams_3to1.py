#make beam table from separate input files (1 per band)
import os
dataDir='/data/Herschel/Calibration/NewBeams/'

beamsOut=TableDataset(description='SPIRE radial beam profiles v3 (Nov 2014)')
beamsOut.addColumn('radius',Column(Float1d.range(1400),description='radius from beam centre (arcsec)'))
for band in ['PSW','PMW','PLW']:
    beamProfIn=asciiTableReader(os.path.join(dataDir,'%s_MCore_9.csv'%band))['c0'].data
    #normalise to 1
    #beamProfIn=beamProfIn/beamProfIn[0]
    beamsOut.addColumn(band,Column(beamProfIn))
    
    
