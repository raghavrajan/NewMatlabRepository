function [] = SimulateHVC_IntroNoteResponses()

% Types of HVC-X projecting neurons based on data
% Equal for all INs - EQUAL
% All but last IN   - ALL BUT LAST
% Last IN           - LAST IN
% Increase firing   - RAMP UP
% Decrease firing   - RAMP DOWN

% I assume that each interneuron gets input from about 