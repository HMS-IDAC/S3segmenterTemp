function S3segmenterWrapperDocker(parametersFile)

javaaddpath('matlabDependencies/bfmatlab/bioformats_package.jar')
javaaddpath('matlabDependencies/yaml/java/snakeyaml-1.9.jar')

[params] = YAML.read(parametersFile);


%Paths that should not change
params.dockerParams.parentPath = '/input';  %docker input folder, mapped by user
params.dockerParams.outputPath = '/output';   %docker output folder, mapped by user
params.dockerParams.modelPath = '/S3segmenter/matlabDependencies';   %.mat, static
params.dockerParams.modelCatPath = '/S3segmenter/matlabDependencies'; %.mat, static

batchS3segmenterWrapper(params.dockerParams.parentPath,params)