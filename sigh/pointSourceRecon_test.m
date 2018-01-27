%% Script to compare axial and full recon (no attenuation correction)
path = '~/Dropbox (RefleXion Medical)/AlphaOne/mockdb-petExp-oct-2017/';
D = rxmNewDeliverySession(path, true);
D.sim.version = 'rxone';
deliveryDeviceId = 'reflexion-alpha-pointsource';
D = setDeliverySession(D, deliveryDeviceId);
D = setPetNormalization( D );

% filename = [path '/' 'cases' '/' 'ImageSeries' '/' '223312-5951969424c477c44cd33216_prescanLor.coincmsg_rxone'];
filename = [path '/' 'regressionTest' '/' '115454-59e55668c7d67532b3fa25b3_prescanLor.coincmsg_rxone'];
D.prescan.lors = readPetLOR(filename,'coincmsg_rxone');
D.prescan.lors = addGeomInfoToLOR( D.prescan.lors );
D.prescan.lors = addNormCorrToLOR(D.prescan.lors, D.norm, D.sinoParam);

%% perform PET recon w/ full mode
protocol = struct(...
    'xyzPixelSzDcmMm', [1 1 2] ...
    , 'imagingExtentDcmMm', [-30 30] ...
    , 'xyzOriginDcmMm', [-255.5 -255.5 -30] ...
    , 'mode', 'full' ...
    , 'interpTyp', 'pchip' ...
    , 'smoothing', 0 ...
    );
orientation = 'PRONE';
eq_rel = 'FEET-FIRST';
dcm_isocenter_xyz_mm = [0 0 0];
treatmentPos = treatmentPosition( orientation, eq_rel, dcm_isocenter_xyz_mm );
mtx = treatmentPos.patientFrameOfReferenceToEquipmentMappingMatrix;
reconParam = setReconParam(D, protocol, mtx);
assert((single(abs(D.norm.normAxialRebinGrid.spacing)) - single(abs(reconParam.imgGrid.xyzSpacing(2))))==0);  %cannot proceed ahead as the rebin eff was generated on wrong grid
D.prescan.reconParam = reconParam;
setupPositionMm = [0 83 0];
D.prescan.lors = discretizeTablePos(D.prescan.lors, D.prescan.reconParam.imgGrid, setupPositionMm);
D.acf = ones(D.sinoParam.nbins, D.sinoParam.nangles, D.prescan.reconParam.ny);
D.prescan.reconParam.RC = 0;
[D.prescan.pet, ~, D.prescan.sino] = runPETRecon(D.prescan.lors, [], D.sinoParam, D.prescan.reconParam, D.acf);
stats1 = struct(...
    'sumSino', isum(D.prescan.sino) ...
    , 'sumImg', isum(D.prescan.pet.Vol) ...
    , 'meanBkgd', imean(D.prescan.pet.Vol(56-10:56+10,46-10:46+10,8-2:8+2)) ...
    , 'stdBkgd', std2(D.prescan.pet.Vol(56-10:56+10,46-10:46+10,8-2:8+2)) ...
    , 'cov', std2(D.prescan.pet.Vol(56-10:56+10,46-10:46+10,8-2:8+2)) / imean(D.prescan.pet.Vol(56-10:56+10,46-10:46+10,8-2:8+2)) ...
    );
toGrid = petCoordTransform(D.prescan.pet.imgGrid, D.commissioning, treatmentPos);
figure;ctshow3(toGrid, D.prescan.pet.Vol, [0 -99 2]);colormap(flipud(gray));

%% Add a unit test for measuring resolution
[PetImg FWHM50 FWHM10] = measure_fwhm(D.prescan.pet);

asdfsdf
