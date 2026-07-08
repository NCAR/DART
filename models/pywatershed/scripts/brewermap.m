function [map,num,typ,scheme] = brewermap(N,scheme) %#ok<*ISMAT>
% The complete selection of ColorBrewer colorschemes/palettes (RGB colormaps).
%
% (c) 2014-2022 Stephen Cobeldick
%
% Returns any RGB colormap from the ColorBrewer colorschemes, especially
% intended for mapping and plots with attractive, distinguishable colors.
%
%%% Basic Syntax:
% brewermap() % print summary
% map = brewermap(N,scheme)
%
%%% Preset Syntax:
% old = brewermap(scheme)
% map = brewermap()
% map = brewermap(N)
%
% [...,num,typ] = brewermap(...)
%
%% Colorschemes %%
%
% This product includes color specifications and designs developed by Cynthia Brewer.
% See the ColorBrewer website for further information about each colorscheme,
% colour-blind suitability, licensing, and citations: http://colorbrewer.org/
% Each colorscheme is defined by a set of hand-picked RGB values (nodes).
%
% To reverse the colormap sequence prefix the colorscheme name with '-'.
%
% Diverging | Qualitative | Sequential
% ----------|-------------|------------------
%  BrBG     |  Accent     |  Blues    PuBuGn
%  PiYG     |  Dark2      |  BuGn     PuRd
%  PRGn     |  Paired     |  BuPu     Purples
%  PuOr     |  Pastel1    |  GnBu     RdPu
%  RdBu     |  Pastel2    |  Greens   Reds
%  RdGy     |  Set1       |  Greys    YlGn
%  RdYlBu   |  Set2       |  OrRd     YlGnBu
%  RdYlGn   |  Set3       |  Oranges  YlOrBr
%  Spectral |             |  PuBu     YlOrRd
%
% If <N> is greater than the requested colorscheme's defining nodes then:
%  - Diverging and Sequential colorschemes are interpolated in Lab colorspace.
%  - Qualitative colorschemes repeat the nodes (i.e. just like LINES does).
% Else:
%  - Exact values from the ColorBrewer colorschemes are returned.
%
%% Examples %%
%
%%% New colors for the COLORMAP example:
% >> S = load('spine');
% >> image(S.X)
% >> colormap(brewermap([],"YlGnBu"))
%
%%% New colors for the SURF example:
% >> [X,Y,Z] = peaks(30);
% >> surfc(X,Y,Z)
% >> colormap(brewermap([],'RdYlGn'))
% >> axis([-3,3,-3,3,-10,5])
%
%%% Plot a colorscheme's RGB values:
% >> rgbplot(brewermap(NaN, 'Blues')) % standard
% >> rgbplot(brewermap(NaN,'-Blues')) % reversed
%
%%% View information about a colorscheme:
% >> [~,num,typ] = brewermap(NaN,'Paired')
% num = 12
% typ = 'Qualitative'
%
%%% Multi-line plot using matrices:
% >> N = 6;
% >> axes('ColorOrder',brewermap(N,'Pastel2'),'NextPlot','replacechildren')
% >> X = linspace(0,pi*3,1000);
% >> Y = bsxfun(@(x,n)n*sin(x+2*n*pi/N), X(:), 1:N);
% >> plot(X,Y, 'linewidth',4)
%
%%% Multi-line plot in a loop:
% set(0,'DefaultAxesColorOrder',brewermap(NaN,'Accent'))
% N = 6;
% X = linspace(0,pi*3,1000);
% Y = bsxfun(@(x,n)n*sin(x+2*n*pi/N), X(:), 1:N);
% for n = 1:N
%     plot(X(:),Y(:,n), 'linewidth',4);
%     hold all
% end
%
%% Input and Output Arguments %%
%
%%% Inputs:
% N = NumericScalar, N>=0, an integer to specify the colormap length.
%   =  [], same length as the current figure's colormap (see COLORMAP).
%   = NaN, same length as the defining RGB nodes (useful for line ColorOrder).
% scheme = CharRowVector or StringScalar, a ColorBrewer colorscheme name.
%
%%% Outputs:
% map = NumericMatrix, size Nx3, a colormap of RGB values between 0 and 1.
% num = NumericVector, the number of nodes defining the ColorBrewer colorscheme.
% typ = CharRowVector, the colorscheme type: 'Diverging'/'Qualitative'/'Sequential'.
%
% See also BREWERMAP_PLOT BREWERMAP_VIEW PRESET_COLORMAP CUBEHELIX MAXDISTCOLOR
% LBMAP PARULA LINES RGBPLOT COLORMAP COLORBAR PLOT PLOT3 AXES SET CONTOURF

%% Input Wrangling %%
%
persistent bmc scm txt
%
if isempty(bmc)
	bmc = bmColors();
end
%
if nargin==0
	N = [];
end
%
err = 'First input <N> must be a real positive scalar numeric or [] or NaN.';
%
if nargout==0 && nargin==0
	hdr = {   'Type'; 'Scheme'; 'Nodes'};
	tsn = [{bmc.typ};{bmc.str};{bmc.num}];
	fprintf('%-12s %-9s %s\n',hdr{:});
	fprintf('%-12s %-9s %u\n',tsn{:});
	return
elseif isnumeric(N)
	if isequal(N,[])
		% Default N is the same as MATLAB colormaps:
		N = cmDefaultN();
	else
		assert(isscalar(N),...
			'SC:brewermap:N:NotScalarNumeric',err)
		assert(isnan(N) || isreal(N) && isfinite(N) && fix(N)==N && N>=0,...
			'SC:brewermap:N:NotWholeRealNotNaN',err)
	end
	if nargin<2
		assert(~isempty(scm),...
			'SC:colorbrewer:scheme:NotPreset',...
			'Colorscheme must be preset before this call: BREWERMAP(SCHEME)')
		scheme = scm;
	else
		scheme = bm1s2c(scheme);
		assert(ischar(scheme) && ndims(scheme)==2 && size(scheme,1)==1,...
			'SC:brewermap:scheme:NotText',...
			'Input <scheme> must be a character vector or string scalar.')
	end
else % preset
	scheme = bm1s2c(N);
	assert(ischar(scheme) && ndims(scheme)==2 && size(scheme,1)==1,...
		'SC:brewermap:N:NotText',...
		'To preset the scheme <N> must be a character vector or string scalar.')
	if strcmpi(scheme,'list')
		map = {bmc.str};
		num = [bmc.num];
		typ = {bmc.typ};
		return
	end
end
%
isr = strncmp(scheme,'-',1) || strncmp(scheme,'*',1);
isd = strncmp(scheme,'+',1) || isr; % direction
ids = strcmpi(scheme(1+isd:end),{bmc.str});
assert(any(ids),...
	'SC:brewermap:scheme:UnknownColorscheme',...
	'Unknown colorscheme name: "%s"',scheme)
%
num = bmc(ids).num;
typ = bmc(ids).typ;
%
if ~isnumeric(N) % preset
	map = txt;
	txt = N;
	scm = scheme;
	return
elseif N==0
	map = nan(0,3);
	return
elseif isnan(N)
	N = num;
end
%
% Downsample:
[idx,itp] = bmIndex(N,num,typ);
%
map = bmc(ids).rgb(idx,:)/255;
%
% Interpolate:
if itp
	%
	M = [... High-precision sRGB to XYZ matrix:
		0.4124564,0.3575761,0.1804375;...
		0.2126729,0.7151522,0.0721750;...
		0.0193339,0.1191920,0.9503041];
	% Source: http://brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
	%
	wpt = [0.95047,1,1.08883]; % D65
	%
	map = bmRGB2Lab(map,M,wpt); % optional
	%
	% Extrapolate a small amount beyond the end nodes:
	%ido = linspace(0,num+1,N+2);
	%ido = ido(2:end-1);
	% Interpolation completely within the end nodes:
	ido = linspace(1,num,N);
	%
	switch typ
		case 'Diverging'
			mid = ceil(num/2);
			ida =   1:mid;
			idz = mid:num;
			map = [...
				interp1(ida,map(ida,:),ido(ido<=mid),'pchip');...
				interp1(idz,map(idz,:),ido(ido>mid),'pchip')];
		case 'Sequential'
			map = interp1(1:num,map,ido,'pchip');
		otherwise
			error('SC:brewermap:NoInterp','Cannot interpolate this type.')
	end
	%
	map = bmLab2RGB(map,M,wpt); % optional
	%
end
%
% Limit output range:
map = max(0,min(1,map));
%
% Reverse row order:
if isr
	map = map(end:-1:1,:);
end
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%brewermap
function N = cmDefaultN()
% Get the colormap size from the current figure or default colormap.
try
	F = get(groot,'CurrentFigure');
catch %#ok<CTCH> pre HG2
	N = size(get(gcf,'colormap'),1);
	return
end
if isempty(F)
	N = size(get(groot,'DefaultFigureColormap'),1);
else
	N = size(F.Colormap,1);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%cmDefaultN
function arr = bm1s2c(arr)
% If scalar string then extract the character vector, otherwise data is unchanged.
if isa(arr,'string') && isscalar(arr)
	arr = arr{1};
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%bm1s2c
function lab = bmRGB2Lab(rgb,M,wpt)
% Convert a matrix of sRGB values to Lab.
%applycform(rgb,makecform('srgb2lab','AdaptedWhitePoint',wpt))
% RGB2XYZ:
xyz = bmGammaInv(rgb) * M.';
% XYZ2Lab:
xyz = bsxfun(@rdivide,xyz,wpt);
idx = xyz>(6/29)^3;
F = idx.*(xyz.^(1/3)) + ~idx.*(xyz*(29/6)^2/3+4/29);
lab(:,2:3) = bsxfun(@times,[500,200],F(:,1:2)-F(:,2:3));
lab(:,1) = 116*F(:,2) - 16;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%bmRGB2Lab
function rgb = bmGammaInv(rgb)
% Inverse gamma correction of sRGB data.
idx = rgb <= 0.04045;
rgb(idx) = rgb(idx) / 12.92;
rgb(~idx) = real(((rgb(~idx) + 0.055) / 1.055).^2.4);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%bmGammaInv
function rgb = bmLab2RGB(lab,M,wpt)
% Convert a matrix of Lab values to sRGB.
%applycform(lab,makecform('lab2srgb','AdaptedWhitePoint',wpt))
% Lab2XYZ
tmp = bsxfun(@rdivide,lab(:,[2,1,3]),[500,Inf,-200]);
tmp = bsxfun(@plus,tmp,(lab(:,1)+16)/116);
idx = tmp>(6/29);
tmp = idx.*(tmp.^3) + ~idx.*(3*(6/29)^2*(tmp-4/29));
xyz = bsxfun(@times,tmp,wpt);
% XYZ2RGB
rgb = max(0,min(1, bmGammaCor(xyz / M.')));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%cbLab2RGB
function rgb = bmGammaCor(rgb)
% Gamma correction of sRGB data.
idx = rgb <= 0.0031308;
rgb(idx) = 12.92 * rgb(idx);
rgb(~idx) = real(1.055 * rgb(~idx).^(1/2.4) - 0.055);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%bmGammaCor
function [idx,itp] = bmIndex(N,num,typ)
% Ensure exactly the same colors as the online ColorBrewer colorschemes.
%
itp = N>num;
switch typ
	case 'Qualitative'
		itp = false;
		idx = 1+mod(0:N-1,num);
	case 'Diverging'
		switch N
			case 1 % extrapolated
				idx = 8;
			case 2 % extrapolated
				idx = [4,12];
			case 3
				idx = [5,8,11];
			case 4
				idx = [3,6,10,13];
			case 5
				idx = [3,6,8,10,13];
			case 6
				idx = [2,5,7,9,11,14];
			case 7
				idx = [2,5,7,8,9,11,14];
			case 8
				idx = [2,4,6,7,9,10,12,14];
			case 9
				idx = [2,4,6,7,8,9,10,12,14];
			case 10
				idx = [1,2,4,6,7,9,10,12,14,15];
			otherwise
				idx = [1,2,4,6,7,8,9,10,12,14,15];
		end
	case 'Sequential'
		switch N
			case 1 % extrapolated
				idx = 6;
			case 2 % extrapolated
				idx = [4,8];
			case 3
				idx = [3,6,9];
			case 4
				idx = [2,5,7,10];
			case 5
				idx = [2,5,7,9,11];
			case 6
				idx = [2,4,6,7,9,11];
			case 7
				idx = [2,4,6,7,8,10,12];
			case 8
				idx = [1,3,4,6,7,8,10,12];
			otherwise
				idx = [1,3,4,6,7,8,10,11,13];
		end
	otherwise
		error('SC:brewermap:UnknownType','Unknown type string.')
end
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%bmIndex
function bmc = bmColors()
% Return a structure of all colorschemes: name, scheme type, RGB values, number of nodes.
% Order: first sort by <typ>, then case-insensitive sort by <str>:
bmc(35).str = 'YlOrRd';
bmc(35).typ = 'Sequential';
bmc(35).rgb = [255,255,204;255,255,178;255,237,160;254,217,118;254,204,92;254,178,76;253,141,60;252,78,42;240,59,32;227,26,28;189,0,38;177,0,38;128,0,38];
bmc(34).str = 'YlOrBr';
bmc(34).typ = 'Sequential';
bmc(34).rgb = [255,255,229;255,255,212;255,247,188;254,227,145;254,217,142;254,196,79;254,153,41;236,112,20;217,95,14;204,76,2;153,52,4;140,45,4;102,37,6];
bmc(33).str = 'YlGnBu';
bmc(33).typ = 'Sequential';
bmc(33).rgb = [255,255,217;255,255,204;237,248,177;199,233,180;161,218,180;127,205,187;65,182,196;29,145,192;44,127,184;34,94,168;37,52,148;12,44,132;8,29,88];
bmc(32).str = 'YlGn';
bmc(32).typ = 'Sequential';
bmc(32).rgb = [255,255,229;255,255,204;247,252,185;217,240,163;194,230,153;173,221,142;120,198,121;65,171,93;49,163,84;35,132,67;0,104,55;0,90,50;0,69,41];
bmc(31).str = 'Reds';
bmc(31).typ = 'Sequential';
bmc(31).rgb = [255,245,240;254,229,217;254,224,210;252,187,161;252,174,145;252,146,114;251,106,74;239,59,44;222,45,38;203,24,29;165,15,21;153,0,13;103,0,13];
bmc(30).str = 'RdPu';
bmc(30).typ = 'Sequential';
bmc(30).rgb = [255,247,243;254,235,226;253,224,221;252,197,192;251,180,185;250,159,181;247,104,161;221,52,151;197,27,138;174,1,126;122,1,119;122,1,119;73,0,106];
bmc(29).str = 'Purples';
bmc(29).typ = 'Sequential';
bmc(29).rgb = [252,251,253;242,240,247;239,237,245;218,218,235;203,201,226;188,189,220;158,154,200;128,125,186;117,107,177;106,81,163;84,39,143;74,20,134;63,0,125];
bmc(28).str = 'PuRd';
bmc(28).typ = 'Sequential';
bmc(28).rgb = [247,244,249;241,238,246;231,225,239;212,185,218;215,181,216;201,148,199;223,101,176;231,41,138;221,28,119;206,18,86;152,0,67;145,0,63;103,0,31];
bmc(27).str = 'PuBuGn';
bmc(27).typ = 'Sequential';
bmc(27).rgb = [255,247,251;246,239,247;236,226,240;208,209,230;189,201,225;166,189,219;103,169,207;54,144,192;28,144,153;2,129,138;1,108,89;1,100,80;1,70,54];
bmc(26).str = 'PuBu';
bmc(26).typ = 'Sequential';
bmc(26).rgb = [255,247,251;241,238,246;236,231,242;208,209,230;189,201,225;166,189,219;116,169,207;54,144,192;43,140,190;5,112,176;4,90,141;3,78,123;2,56,88];
bmc(25).str = 'Oranges';
bmc(25).typ = 'Sequential';
bmc(25).rgb = [255,245,235;254,237,222;254,230,206;253,208,162;253,190,133;253,174,107;253,141,60;241,105,19;230,85,13;217,72,1;166,54,3;140,45,4;127,39,4];
bmc(24).str = 'OrRd';
bmc(24).typ = 'Sequential';
bmc(24).rgb = [255,247,236;254,240,217;254,232,200;253,212,158;253,204,138;253,187,132;252,141,89;239,101,72;227,74,51;215,48,31;179,0,0;153,0,0;127,0,0];
bmc(23).str = 'Greys';
bmc(23).typ = 'Sequential';
bmc(23).rgb = [255,255,255;247,247,247;240,240,240;217,217,217;204,204,204;189,189,189;150,150,150;115,115,115;99,99,99;82,82,82;37,37,37;37,37,37;0,0,0];
bmc(22).str = 'Greens';
bmc(22).typ = 'Sequential';
bmc(22).rgb = [247,252,245;237,248,233;229,245,224;199,233,192;186,228,179;161,217,155;116,196,118;65,171,93;49,163,84;35,139,69;0,109,44;0,90,50;0,68,27];
bmc(21).str = 'GnBu';
bmc(21).typ = 'Sequential';
bmc(21).rgb = [247,252,240;240,249,232;224,243,219;204,235,197;186,228,188;168,221,181;123,204,196;78,179,211;67,162,202;43,140,190;8,104,172;8,88,158;8,64,129];
bmc(20).str = 'BuPu';
bmc(20).typ = 'Sequential';
bmc(20).rgb = [247,252,253;237,248,251;224,236,244;191,211,230;179,205,227;158,188,218;140,150,198;140,107,177;136,86,167;136,65,157;129,15,124;110,1,107;77,0,75];
bmc(19).str = 'BuGn';
bmc(19).typ = 'Sequential';
bmc(19).rgb = [247,252,253;237,248,251;229,245,249;204,236,230;178,226,226;153,216,201;102,194,164;65,174,118;44,162,95;35,139,69;0,109,44;0,88,36;0,68,27];
bmc(18).str = 'Blues';
bmc(18).typ = 'Sequential';
bmc(18).rgb = [247,251,255;239,243,255;222,235,247;198,219,239;189,215,231;158,202,225;107,174,214;66,146,198;49,130,189;33,113,181;8,81,156;8,69,148;8,48,107];
bmc(17).str = 'Set3';
bmc(17).typ = 'Qualitative';
bmc(17).rgb = [141,211,199;255,255,179;190,186,218;251,128,114;128,177,211;253,180,98;179,222,105;252,205,229;217,217,217;188,128,189;204,235,197;255,237,111];
bmc(16).str = 'Set2';
bmc(16).typ = 'Qualitative';
bmc(16).rgb = [102,194,165;252,141,98;141,160,203;231,138,195;166,216,84;255,217,47;229,196,148;179,179,179];
bmc(15).str = 'Set1';
bmc(15).typ = 'Qualitative';
bmc(15).rgb = [228,26,28;55,126,184;77,175,74;152,78,163;255,127,0;255,255,51;166,86,40;247,129,191;153,153,153];
bmc(14).str = 'Pastel2';
bmc(14).typ = 'Qualitative';
bmc(14).rgb = [179,226,205;253,205,172;203,213,232;244,202,228;230,245,201;255,242,174;241,226,204;204,204,204];
bmc(13).str = 'Pastel1';
bmc(13).typ = 'Qualitative';
bmc(13).rgb = [251,180,174;179,205,227;204,235,197;222,203,228;254,217,166;255,255,204;229,216,189;253,218,236;242,242,242];
bmc(12).str = 'Paired';
bmc(12).typ = 'Qualitative';
bmc(12).rgb = [166,206,227;31,120,180;178,223,138;51,160,44;251,154,153;227,26,28;253,191,111;255,127,0;202,178,214;106,61,154;255,255,153;177,89,40];
bmc(11).str = 'Dark2';
bmc(11).typ = 'Qualitative';
bmc(11).rgb = [27,158,119;217,95,2;117,112,179;231,41,138;102,166,30;230,171,2;166,118,29;102,102,102];
bmc(10).str = 'Accent';
bmc(10).typ = 'Qualitative';
bmc(10).rgb = [127,201,127;190,174,212;253,192,134;255,255,153;56,108,176;240,2,127;191,91,23;102,102,102];
bmc(09).str = 'Spectral';
bmc(09).typ = 'Diverging';
bmc(09).rgb = [158,1,66;213,62,79;215,25,28;244,109,67;252,141,89;253,174,97;254,224,139;255,255,191;230,245,152;171,221,164;153,213,148;102,194,165;43,131,186;50,136,189;94,79,162];
bmc(08).str = 'RdYlGn';
bmc(08).typ = 'Diverging';
bmc(08).rgb = [165,0,38;215,48,39;215,25,28;244,109,67;252,141,89;253,174,97;254,224,139;255,255,191;217,239,139;166,217,106;145,207,96;102,189,99;26,150,65;26,152,80;0,104,55];
bmc(07).str = 'RdYlBu';
bmc(07).typ = 'Diverging';
bmc(07).rgb = [165,0,38;215,48,39;215,25,28;244,109,67;252,141,89;253,174,97;254,224,144;255,255,191;224,243,248;171,217,233;145,191,219;116,173,209;44,123,182;69,117,180;49,54,149];
bmc(06).str = 'RdGy';
bmc(06).typ = 'Diverging';
bmc(06).rgb = [103,0,31;178,24,43;202,0,32;214,96,77;239,138,98;244,165,130;253,219,199;255,255,255;224,224,224;186,186,186;153,153,153;135,135,135;64,64,64;77,77,77;26,26,26];
bmc(05).str = 'RdBu';
bmc(05).typ = 'Diverging';
bmc(05).rgb = [103,0,31;178,24,43;202,0,32;214,96,77;239,138,98;244,165,130;253,219,199;247,247,247;209,229,240;146,197,222;103,169,207;67,147,195;5,113,176;33,102,172;5,48,97];
bmc(04).str = 'PuOr';
bmc(04).typ = 'Diverging';
bmc(04).rgb = [127,59,8;179,88,6;230,97,1;224,130,20;241,163,64;253,184,99;254,224,182;247,247,247;216,218,235;178,171,210;153,142,195;128,115,172;94,60,153;84,39,136;45,0,75];
bmc(03).str = 'PRGn';
bmc(03).typ = 'Diverging';
bmc(03).rgb = [64,0,75;118,42,131;123,50,148;153,112,171;175,141,195;194,165,207;231,212,232;247,247,247;217,240,211;166,219,160;127,191,123;90,174,97;0,136,55;27,120,55;0,68,27];
bmc(02).str = 'PiYG';
bmc(02).typ = 'Diverging';
bmc(02).rgb = [142,1,82;197,27,125;208,28,139;222,119,174;233,163,201;241,182,218;253,224,239;247,247,247;230,245,208;184,225,134;161,215,106;127,188,65;77,172,38;77,146,33;39,100,25];
bmc(01).str = 'BrBG';
bmc(01).typ = 'Diverging';
bmc(01).rgb = [84,48,5;140,81,10;166,97,26;191,129,45;216,179,101;223,194,125;246,232,195;245,245,245;199,234,229;128,205,193;90,180,172;53,151,143;1,133,113;1,102,94;0,60,48];
% number of nodes:
for k = 1:numel(bmc)
	switch bmc(k).typ
		case 'Diverging'
			bmc(k).num = 11;
		case 'Qualitative'
			bmc(k).num = size(bmc(k).rgb,1);
		case 'Sequential'
			bmc(k).num = 9;
		otherwise
			error('SC:brewermap:UnknownType','Unknown type string.')
	end
end
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%bmColors
%
% Code and Implementation:
% Copyright (c) 2014-2022 Stephen Cobeldick
% Color Values Only:
% Copyright (c) 2002 Cynthia Brewer, Mark Harrower, and The Pennsylvania State University.
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
% http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and limitations under the License.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions as source code must retain the above copyright notice, this
% list of conditions and the following disclaimer.
%
% 2. The end-user documentation included with the redistribution, if any, must
% include the following acknowledgment: "This product includes color
% specifications and designs developed by Cynthia Brewer
% (http://colorbrewer.org/)." Alternately, this acknowledgment may appear in the
% software itself, if and wherever such third-party acknowledgments normally appear.
%
% 4. The name "ColorBrewer" must not be used to endorse or promote products
% derived from this software without prior written permission. For written
% permission, please contact Cynthia Brewer at cbrewer@psu.edu.
%
% 5. Products derived from this software may not be called "ColorBrewer", nor
% may "ColorBrewer" appear in their name, without prior written permission of Cynthia Brewer.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%license
