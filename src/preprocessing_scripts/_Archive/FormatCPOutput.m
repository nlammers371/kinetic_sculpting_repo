addpath('./utilities')
%------------------------Import Compiled Particles------------------------%
project_list = {'2017-02-08-eve23_4_noCur',...
                '2017-05-13-A200umP_eve_5uW'...
                '2017-05-14-A150umP_eve_5uW'}; 
%project_list = {'2017-05-13-A200umP_eve_5uW'};
datapath = 'D:\Data\Nick\LivemRNA\LivemRNAFISH\Data\DynamicsResults\';
outpath = 'D:\Data\Nick\AugustoData\ForMike';
if exist(outpath) ~= 7
    mkdir(outpath)
end

for file = project_list
   
    traces = load([datapath file{:} '\CompiledParticles.mat']);
    particles = traces.CompiledParticles;
    particles = particles([particles.nc]==14);
    time = traces.ElapsedTime*60 - traces.ElapsedTime(traces.nc14);
    for j = 1:length(particles)
        pFrames = particles(j).Frame;
        FrameIndex = min(pFrames):max(pFrames);
        particles(j).Time = time(FrameIndex);
    end
    raw_traces = particles;
    nucleus_data = load([datapath file{:} '\' file{:} '_lin.mat']);

    n_particles = length(raw_traces);
    %Get nuclear positions
    i_exclude = [];
    for i = 1:n_particles
        n_id = raw_traces(i).schnitz;
        xPos = nucleus_data.schnitzcells(n_id).cenx;
        yPos = nucleus_data.schnitzcells(n_id).ceny;
        nuc_Frames = nucleus_data.schnitzcells(n_id).frames;
        if length(raw_traces(i).xPos) ~= length(xPos(ismember(nuc_Frames,raw_traces(i).Frame)))
            i_exclude = [i_exclude i];     
        else
            xRel = raw_traces(i).xPos - xPos(ismember(nuc_Frames,raw_traces(i).Frame));
            yRel = raw_traces(i).yPos - yPos(ismember(nuc_Frames,raw_traces(i).Frame));
            raw_traces(i).xRel = xRel;
            raw_traces(i).yRel = yRel;
        end
    end
    ind = 1:length(raw_traces);
    raw_traces = raw_traces(~ismember(ind,i_exclude));
    traces_out = struct;
    for i = 1:length(raw_traces)
        particle = raw_traces(i);
        Fluo = particle.Fluo;
        pFrames = particle.Frame;

        FrameIndex = pFrames - min(pFrames) + 1;
        FramesFull = min(pFrames):max(pFrames);

        FluoAugmented = zeros(1,length(FramesFull));
        FluoAugmented(FrameIndex) = Fluo(:);
        %Zero any negative fluorecence values (16 points total)
        FluoAugmented(FluoAugmented < 0) = 0;
      
        %Interpolate AP and nc rel pos info
        ap_orig = raw_traces(i).APPos;
        ap_new = interp1(pFrames, ap_orig, FramesFull);
        xPosNew = interp1(pFrames,raw_traces(i).xRel,FramesFull);
        yPosNew = interp1(pFrames,raw_traces(i).yRel,FramesFull);
        %Store Updated info
        traces_out(i).AP = ap_new;
        traces_out(i).xRel = xPosNew;
        traces_out(i).yRel = yPosNew;
        traces_out(i).Fluo = FluoAugmented;
        traces_out(i).Time = raw_traces(i).Time;
        traces_out(i).OriginalParticle = linspace(raw_traces(i).OriginalParticle,...
            raw_traces(i).OriginalParticle,length(raw_traces(i).Time));
    end

    output = [[traces_out.OriginalParticle]', [traces_out.Time]',...
              [traces_out.Fluo]', ...
              [traces_out.AP]', [traces_out.xRel]', [traces_out.yRel]'];
    
    header = {'ParticleID','Seconds','Fluo','AP','xRel','yRel'};
    csvwrite_with_headers([outpath '\' file{:} '_longform.csv'], ...
                       output, header);      
   
end