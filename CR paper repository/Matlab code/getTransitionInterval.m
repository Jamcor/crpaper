function [TransitionIntervals,StartTime]=getTransitionInterval(listmode,ndx,event_type)
% gets transition interval for the event type 
% measures time from the start of the experiment rather than cell age
% listmode is the wavelength field
% clones is a vector of clone IDs
% is a vector of track IDs (same length as clone IDs)
switch(event_type)
    case 'GFP_on'
        [TransitionIntervals,StartTime]=getGFPTime(listmode,ndx);
    case 'Death'
        [TransitionIntervals,StartTime]=getRightCensoredTime(listmode,ndx);
    case 'Right_Censored'
        [TransitionIntervals,StartTime]=getRightCensoredTime(listmode,ndx);
    case 'Division'
        [TransitionIntervals,StartTime]=getRightCensoredTime(listmode,ndx);
    otherwise 
        error('Event_type not recognised');
end
 
function [t,StartTime]=getGFPTime(listmode,ndx)
clones=listmode(1,ndx);
t=listmode(10,ndx)+listmode(7,ndx); %add birth time of the cell to the age of transition (gives duration from beginning of experiment until transition) 
StartTime=zeros(size(t));
%first row is clone ID, 2nd row is track number, row 10 is birth time
for i=1:length(clones)
    anc_ndx=(listmode(2,:)==1)&(listmode(1,:)==clones(i)); %get index for the first track for first clone
    anc_birthtime=listmode(10,anc_ndx); %get the birtime of founding cell of clone
    t(i)=t(i)-anc_birthtime;  %will give the time from the birth of founding cell until GFP switches on in progeny (
    StartTime(i)=anc_birthtime; 
    
    %t gives duration from beginning of experiment until transition) 
end
 
function [t,StartTime]=getRightCensoredTime(listmode,ndx)
clones=listmode(1,ndx);
t=listmode(10,ndx)+listmode(8,ndx);
StartTime=zeros(size(t));
%first row is clone ID, 2nd row is birth time
for i=1:length(clones)
    anc_ndx=(listmode(2,:)==1)&(listmode(1,:)==clones(i));
    anc_birthtime=listmode(10,anc_ndx);
    t(i)=t(i)-anc_birthtime; 
    StartTime(i)=anc_birthtime;
end
 