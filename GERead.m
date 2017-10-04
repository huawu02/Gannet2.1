function [ MRS_struct ] = GERead(MRS_struct, fname)
ii=MRS_struct.ii;
MRS_struct.p.global_rescale=1;
%121106 RAEE moving code from GAnnetLoad to GERead (to match other file
%formats and tidy things up some.
%RTN edits to accommodate Noeske version RAEE 141007
%160916: MM & RTN edits to accommodate different encoding schemes
% HW: should use ieee-le to read pfile
fid = fopen(fname,'r', 'ieee-le'); 
if fid == -1
    tmp = [ 'Unable to locate Pfile ' fname ];
    disp(tmp);
    return;
end
% return error message if unable to read file type.
% Determine size of Pfile header based on Rev number
fseek(fid, 0, 'bof');
f_hdr_value = fread(fid, 1, 'real*4');
rdbm_rev_num = f_hdr_value(1);
if( rdbm_rev_num == 7.0 )
    pfile_header_size = 39984;  % LX
elseif ( rdbm_rev_num == 8.0 )
    pfile_header_size = 60464;  % Cardiac / MGD
elseif (( rdbm_rev_num > 5.0 ) && (rdbm_rev_num < 6.0))
    pfile_header_size = 39940;  % Signa 5.5
else
    % In 11.0 and later the header and data are stored as little-endian
    fclose(fid);
    fid = fopen(fname,'r', 'ieee-le');
    fseek(fid, 0, 'bof');
    f_hdr_value = fread(fid, 1, 'real*4');
    if (f_hdr_value == 9.0)  % 11.0 product release
        pfile_header_size= 61464;
    elseif (f_hdr_value == 11.0)  % 12.0 product release
        pfile_header_size= 66072;
    elseif (f_hdr_value > 11.0) && (f_hdr_value < 26.0)  % 14.0 and later, before 26.0
        fseek(fid, 1468, 'bof');
        pfile_header_size = fread(fid,1,'integer*4');
    elseif (f_hdr_value >= 26.0)    % HW: 26.0 and later
        pfile_header_size = 213684;
    else
        sprintf('Invalid Pfile header revision: %f', f_hdr_value );
        return;
    end
end

MRS_struct.p.GE.rdbm_rev_num = floor(f_hdr_value(1)*1000)/1000; % MM (170118)
chkRev = [16, 24, 26]; % GERead mods tested with these revisions only. HW: add rev 26 support
if ~any(floor(MRS_struct.p.GE.rdbm_rev_num) == chkRev)
    warning('!!! GERead not fully functional with header revision number %d !!!', MRS_struct.p.GE.rdbm_rev_num);
end

if MRS_struct.p.GE.rdbm_rev_num == 26.002     % HW: add rev 26 support
    fseek(fid, 140, 'bof');     npasses                                 = fread(fid, 1, 'int16');
    fseek(fid, 144, 'bof');     nslices                                 = fread(fid, 1, 'int16');
    fseek(fid, 146, 'bof');     nechoes                                 = fread(fid, 1, 'int16');
    fseek(fid, 148, 'bof');     nex                                     = fread(fid, 1, 'int16');           % rdb_hdr.navs
    fseek(fid, 150, 'bof');     nframes                                 = fread(fid, 1, 'int16');
    fseek(fid, 158, 'bof');     point_size                              = fread(fid, 1, 'int16');
    fseek(fid, 178, 'bof');     MRS_struct.p.npoints                    = fread(fid, 1, 'int16');           % rdb_hdr.da_xres
    fseek(fid, 180, 'bof');     MRS_struct.p.nrows                      = fread(fid, 1, 'int16');           % rdb_hdr.da_yres
    fseek(fid, 264, 'bof');     start_recv                              = fread(fid, 1, 'int16');           % rdb_hdr.dab[0]
    fseek(fid, 266, 'bof');     stop_recv                               = fread(fid, 1, 'int16');           % rdb_hdr.dab[1]
    fseek(fid, 280, 'bof');     MRS_struct.p.sw                         = fread(fid, 1, 'float32');         % rdb_hdr.user0
    if nechoes ~= 1
        fseek(fid, 296, 'bof');     dataframes                          = fread(fid, 1, 'float32') / nex;   % rdb_hdr.user4
        fseek(fid, 356, 'bof');     refframes                           = fread(fid, 1, 'float32');         % rdb_hdr.user19
    end
    fseek(fid, 488, 'bof');     MRS_struct.p.LarmorFreq                 = fread(fid, 1, 'uint32') / 1e7;    % rdb_hdr.ps_mps_freq
    fseek(fid, 199244, 'bof');  MRS_struct.p.TE(ii)                     = fread(fid, 1, 'int32') / 1e3;     % image_hdr.te
    fseek(fid, 199236, 'bof');  MRS_struct.p.TR(ii)                     = fread(fid, 1, 'int32') / 1e3;
    fseek(fid, 198568, 'bof');  MRS_struct.p.voxdim(ii,:)               = fread(fid, 3, 'float32')';        % image_hdr.user8 : image_hdr.user10
    fseek(fid, 198612, 'bof');  MRS_struct.p.GE.editRF.waveform(ii)     = fread(fid, 1, 'float32');         % image_hdr.user19
    fseek(fid, 198616, 'bof');  MRS_struct.p.GE.editRF.freq_Hz(ii,:)    = fread(fid, 2, 'float32')';        % image_hdr.user20 : 21
    fseek(fid, 198624, 'bof');  MRS_struct.p.GE.editRF.dur(ii)          = fread(fid, 1, 'float32') / 1e3;   % image_hdr.user22
    % get the exam/series number, for naming the output file
    fseek(fid, 199552, 'bof');  MRS_struct.p.ex_no                      = fread(fid, 1, 'uint16');          % image_hdr.im_exno
    fseek(fid, 199564, 'bof');  MRS_struct.p.se_no                      = fread(fid, 1, 'int16');           % image_hdr.im_seno

else     % if any(MRS_struct.p.rdbm_rev_num == [16, 24])
    % Read header information
    fseek(fid, 0, 'bof');
    hdr_value = fread(fid, 102, 'integer*2');
    % RTN - read rhuser
    fseek(fid, 0, 'bof');
    f_hdr_value = fread(fid, 74, 'real*4');
     % RTN (170118): Find center frequency
    fseek(fid, 0, 'bof');
    i_hdr_value = fread(fid, 102+9, 'integer*4');
    MRS_struct.p.LarmorFreq = i_hdr_value(102+5)/1e7;
    MRS_struct.p.sw = f_hdr_value(55); % MM (160916)

    npasses = hdr_value(33);
    nslices = hdr_value(35);
    nechoes = hdr_value(36);
    %RTN - number of phase cycles
    nex = hdr_value(37);
    nframes = hdr_value(38);
    point_size = hdr_value(42);
    MRS_struct.p.npoints = hdr_value(52);
    MRS_struct.p.nrows = hdr_value(53);
    % rc_xres = hdr_value(54);
    % rc_yres = hdr_value(55);
    start_recv = hdr_value(101);
    stop_recv = hdr_value(102);
    if nechoes ~= 1
        dataframes = f_hdr_value(59)/nex;
        refframes = f_hdr_value(74);
    end
    % MM (170118): Find TE/TR
    fseek(fid, 1468, 'bof');
    p_hdr_value = fread(fid, 12, 'integer*4'); % byte offsets to start of sub-header structures
    fseek(fid, p_hdr_value(10), 'bof'); % set position to start of rdb_hdr_image
    t_hdr_value = fread(fid, p_hdr_value(1)-p_hdr_value(10), 'integer*4');
    if isequal(MRS_struct.p.GE.rdbm_rev_num, 16)    
        MRS_struct.p.TE(ii) = t_hdr_value(193)/1e3;
        MRS_struct.p.TR(ii) = t_hdr_value(191)/1e3;
    elseif isequal(MRS_struct.p.GE.rdbm_rev_num, 24)
        MRS_struct.p.TE(ii) = t_hdr_value(267)/1e3;
        MRS_struct.p.TR(ii) = t_hdr_value(265)/1e3;
    end

    % MM (170127): Find voxel dimensions and edit pulse parameters
    fseek(fid, p_hdr_value(8), 'bof'); % set position to start of rdb_hdr_exam. HW: actually these fields are in rdb_hdr_image.user*
    o_hdr_value = fread(fid, p_hdr_value(9)-p_hdr_value(8), 'real*4');
    if isequal(MRS_struct.p.GE.rdbm_rev_num, 16)
        MRS_struct.p.voxdim(ii,:) = o_hdr_value(822:824)';
        MRS_struct.p.GE.editRF.waveform(ii) = o_hdr_value(833);
        MRS_struct.p.GE.editRF.freq_Hz(ii,:) = o_hdr_value(834:835)';
        MRS_struct.p.GE.editRF.dur(ii) = o_hdr_value(836)/1e3;
    elseif isequal(MRS_struct.p.GE.rdbm_rev_num, 24)
        MRS_struct.p.voxdim(ii,:) = o_hdr_value(1228:1230)';
        MRS_struct.p.GE.editRF.waveform(ii) = o_hdr_value(1239);
        MRS_struct.p.GE.editRF.freq_Hz(ii,:) = o_hdr_value(1240:1241)';
        MRS_struct.p.GE.editRF.dur(ii) = o_hdr_value(1242)/1e3;
    end
    
    % get the exam/series number, for naming the output file
    fseek(fid, 148712, 'bof');  MRS_struct.p.ex_no                   = fread(fid, 1, 'uint16');          % image_hdr.im_exno
    fseek(fid, 148724, 'bof');  MRS_struct.p.se_no                   = fread(fid, 1, 'int16');           % image_hdr.im_seno

end

MRS_struct.p.GE.nechoes = nechoes;
MRS_struct.p.GE.NEX = nex;
nreceivers = (stop_recv - start_recv) + 1;
MRS_struct.p.GE.editRF.freq_ppm(ii,:) = (MRS_struct.p.GE.editRF.freq_Hz(ii,:) / MRS_struct.p.LarmorFreq) + 4.68;

% Spectro prescan pfiles
if (MRS_struct.p.npoints == 1) && (MRS_struct.p.nrows == 1)
    MRS_struct.p.npoints = 2048;
end

% Determine number of slices in this Pfile:  this does not work for all cases.
slices_in_pass = nslices/npasses;

% Compute size (in bytes) of each frame, echo and slice
data_elements = MRS_struct.p.npoints*2;
frame_size = data_elements*point_size;
echo_size = frame_size*MRS_struct.p.nrows;
slice_size = echo_size*nechoes;
mslice_size = slice_size*slices_in_pass;
my_slice = 1;
my_echo = 1;
my_frame = 1;

FullData=zeros(nreceivers, MRS_struct.p.npoints, (MRS_struct.p.nrows-my_frame+1)*nechoes); %RTN nechoes multiplication;

%Start to read data into Eightchannel structure.
totalframes=(MRS_struct.p.nrows-my_frame+1)*nechoes; % RTN nechoes mulitply;
MRS_struct.p.nrows=totalframes;
data_elements2 = data_elements*totalframes*nreceivers;

%  % Compute offset in bytes to start of frame.
file_offset = pfile_header_size + ((my_frame-1)*frame_size);

fseek(fid, file_offset, 'bof');

% read data: point_size = 2 means 16 bit data, point_size = 4 means EDR )
if (point_size == 2 )
    raw_data = fread(fid, data_elements2, 'integer*2');
else
    raw_data = fread(fid, data_elements2, 'integer*4');
end

fclose(fid);

% 110303 CJE
% calculate Navg from nframes, 8 water frames, 2 phase cycles
% Needs to be specific to single experiment - for frame rejection
% RTN edits to accommodate Noeske version raee 20141007
% MM (160916): Incorporating more edits from RTN to handle dual-echo data
%              acquired with one of four possible encoding schemes:
%              NEX=2/noadd=0, NEX=2/noadd=1, NEX=8/noadd=0, NEX=8/noadd=1
if (nechoes == 1)
    MRS_struct.p.Navg(ii) = (nframes-8)*2;
    MRS_struct.p.Nwateravg = 8; %moved from MRSGABAinstunits RE 110726
    ShapeData = reshape(raw_data,[2 MRS_struct.p.npoints totalframes nreceivers]);
    ZeroData = ShapeData(:,:,1,:);
    WaterData = ShapeData(:,:,2:9,:);
    FullData = ShapeData(:,:,10:end,:);
    
    totalframes = totalframes-9;
    MRS_struct.p.nrows=totalframes;
    
    Frames_for_Water = 8;
else
    MRS_struct.p.Navg(ii) = dataframes*nex*2; % RTN 2016
    
    if ((dataframes+refframes) ~= nframes)
        mult = nex/2.0; % RTN 2016
        multw = nex; % RTN 2016
        MRS_struct.p.GE.noadd = 1;
        dataframes = dataframes*nex;
        refframes = nframes - dataframes; % refframes*nex; 2015
    else
        mult = nex/2.0; % RTN 2016
        multw = 1.0; % RTN 2016
        MRS_struct.p.GE.noadd = 0;
    end
    
    MRS_struct.p.Nwateravg(ii) = refframes*2;
    
    if (totalframes ~= ((dataframes+refframes+1)*2))
        error('# of totalframes not same as (dataframes+refframes+1)*2');
    end
    ShapeData = reshape(raw_data,[2 MRS_struct.p.npoints totalframes nreceivers]);
    ZeroData = ShapeData(:,:,1,:);
    WaterData = zeros([2 MRS_struct.p.npoints refframes*2 nreceivers]);
    for loop = 1:refframes
        WaterData(:,:,2*loop,:)=(-1)^(MRS_struct.p.GE.noadd*(loop-1))*ShapeData(:,:,1+loop,:) * multw; % RTN 2016
        WaterData(:,:,2*loop-1,:)=(-1)^(MRS_struct.p.GE.noadd*(loop-1))*ShapeData(:,:,totalframes/2+1+loop,:) * multw; % RTN 2016
    end
    FullData = zeros([2 MRS_struct.p.npoints dataframes*2 nreceivers]);
    for loop = 1:dataframes
        FullData(:,:,2*loop,:)=(-1)^(MRS_struct.p.GE.noadd*(loop-1))*ShapeData(:,:,1+refframes+loop,:) * mult; % RTN 2016
        FullData(:,:,2*loop-1,:)=(-1)^(MRS_struct.p.GE.noadd*(loop-1))*ShapeData(:,:,totalframes/2+refframes+1+loop,:) * mult; % RTN 2016
    end
    totalframes=totalframes-refframes*2-2;
    MRS_struct.p.nrows=totalframes;
    Frames_for_Water=refframes*2;
end

FullData = FullData.*repmat([1;1i],[1 MRS_struct.p.npoints totalframes nreceivers]);
FullData = squeeze(sum(FullData,1));
FullData = permute(FullData,[3 1 2]);
WaterData = WaterData.*repmat([1;1i],[1 MRS_struct.p.npoints Frames_for_Water nreceivers]);
WaterData = squeeze(sum(WaterData,1));
WaterData = permute(WaterData,[3 1 2]);
% at this point, FullData(rx_channel, point, average)

% % Old code
% firstpoint=conj(WaterData(:,1,:));
% firstpoint=repmat(firstpoint, [1 MRS_struct.p.npoints 1]);
% % here firstpoint(rx_channel,[], average)
% 
% % CJE March 10 - correct phase of each Water avg independently
% WaterData=WaterData.*firstpoint*MRS_struct.p.global_rescale;
% 
% %Multiply the Eightchannel data by the firstpointvector
% % zeroth order phasing of spectra
% % CJE Nov 09: do global rescaling here too
% % don't really need the phasing step here if performing frame-by-frame phasing
% for receiverloop = 1:nreceivers
%     FullData(receiverloop,:) = FullData(receiverloop,:)*firstpoint(receiverloop,1,1)*MRS_struct.p.global_rescale;
%     % WaterData(receiverloop,:) = WaterData(receiverloop,:)*firstpoint(receiverloop,1,1)*MRS_struct.global_rescale;
% end
% 
% % sum over Rx channels
% FullData = squeeze(sum(FullData,1));
% MRS_struct.fids.data =FullData;
% WaterData = squeeze(sum(WaterData,1));
% MRS_struct.fids.data_water=WaterData;

% MM (170505)
firstpoint_water = conj(WaterData(:,1,:));
channels_scale = squeeze(sqrt(sum(firstpoint_water .* conj(firstpoint_water),1)));
channels_scale = repmat(channels_scale, [1 nreceivers MRS_struct.p.npoints]);
channels_scale = permute(channels_scale, [2 3 1]);
firstpoint_water = repmat(firstpoint_water, [1 MRS_struct.p.npoints 1])./channels_scale;

WaterData = WaterData .* firstpoint_water * MRS_struct.p.global_rescale;
WaterData = squeeze(sum(WaterData,1));
MRS_struct.fids.data_water = WaterData;

% Use first point of water data to phase water-suppressed data
firstpoint = mean(firstpoint_water,3);
firstpoint = repmat(firstpoint, [1 1 size(FullData,3)]);

FullData = FullData .* firstpoint * MRS_struct.p.global_rescale;
FullData = squeeze(sum(FullData,1));
MRS_struct.fids.data = FullData;

% % S/N^2 coil-combination method (MM: 170505)
% phi = unwrap(angle(WaterData(:,1,:)));
% phi = exp(-1i.*phi);
% phi_WaterData = repmat(phi, [1 MRS_struct.p.npoints 1]);
% 
% S = max(abs(WaterData(:,1,:)));
% N = std(WaterData(:,end-100:end,:),[],2);
% w = S./N.^2;
% w_WaterData = repmat(w, [1 MRS_struct.p.npoints 1]);
% 
% WaterData = WaterData .* phi_WaterData .* w_WaterData * MRS_struct.p.global_rescale;
% WaterData = squeeze(sum(WaterData,1));
% MRS_struct.fids.data_water = WaterData;
% 
% phi_FullData = repmat(mean(phi,3), [1 MRS_struct.p.npoints totalframes]);
% w_FullData = repmat(mean(w,3), [1 MRS_struct.p.npoints totalframes]);
% 
% FullData = FullData .* phi_FullData .* w_FullData * MRS_struct.p.global_rescale;
% FullData = squeeze(sum(FullData,1));
% MRS_struct.fids.data = FullData;

%%%%%% end of GE specific load

rescale = 1/1e11; % necessary for GE data or numbers blow up
MRS_struct.fids.data = MRS_struct.fids.data * rescale;
MRS_struct.fids.data_water = MRS_struct.fids.data_water * rescale;

end


