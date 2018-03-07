import * as React from 'react'; 
import {Card, CardActions, CardHeader, CardText} from 'material-ui/Card';
import FlatButton from 'material-ui/FlatButton';

import Paper from 'material-ui/Paper';

import SelectedElements from '../components/SelectedElements';
import ACInput from '../components/AutoComplete';
import MSATableViewer from '../components/MSATableViewer';

import axios from 'axios';
import config from '../config';
  
export interface QueryComponentProps { key: number};
export interface QueryComponentState { selectedElements: Array<any>, alignments: any };

class QueryComponent extends React.Component<QueryComponentProps, QueryComponentState> {
    constructor(props) {
        super(props);

    }

    componentWillMount()
    {
        this.setState({selectedElements: []});
    }

    newElementSelected( newElement )
    {
        this.state.selectedElements.push(newElement);
        this.setState({selectedElements: this.state.selectedElements});
    }

    elementClicked( elementText )
    {
        console.log("Clicked on: " + elementText)
    }

    deleteElement( elementText, i )
    {
        var idx = this.state.selectedElements.indexOf(elementText);

        if (idx >= 0)
        {
            this.state.selectedElements.splice(idx, 1);
        }

        this.setState({selectedElements: this.state.selectedElements})
    }

    prepareResults()
    {
        var self = this;

        axios.post(config.getRestAddress() + "/clustalign", {genes: this.state.selectedElements}, config.axiosConfig)
        .then(function (response) {
          console.log(response.data)

          self.setState({alignments: response.data})

        })
        .catch(function (error) {
          console.log(error)
          self.setState({alignments: {}})
        });


    }

    render()
    {

        var alignResults = [];
        if ((this.state.alignments == null) || (this.state.alignments.length == 0))
        {
            alignResults.push(<p key={0}>No Result Available yet</p>)   ;
        } else {

            var alignKeys = Object.keys(this.state.alignments);
            console.log("AlignKeys");
            console.log(alignKeys);
            
            for (var ki=0; ki < alignKeys.length; ++ki)
            {
                var key = alignKeys[ki];
                console.log("KEY");
                console.log(key);
                
                for (var hi = 0; hi < this.state.alignments[key].length; ++hi)
                {
                    var homCluster = this.state.alignments[key][hi];
                    console.log("HOMCLUSTER");
                    console.log(homCluster);
                    
                    var homClusterAlign = {'id': homCluster['homid'], 'align': []}
                    
                    for (var si=0; si < homCluster['msa'].length; ++si)
                    {
                        var seq = homCluster['msa'][si];
                        var name = seq['org'] + " " + seq['seqid'];

                        homClusterAlign.align.push({org: seq['org'], seqid: seq['seqid'], id: name, alignment: seq['alignment'], xrefs: seq['xrefs'] })
                    }
                    console.log("HOMALIGN");
                    console.log(homClusterAlign);
                    alignResults.push( <MSATableViewer key={alignResults.length} alignments={homClusterAlign}/> )
                }

            }
        }

        return (

            <Card>
                <CardHeader
                title="Show Homologies"
                subtitle=""
                />
                <CardText>

                    <div>
                        <ACInput onElementSelected={this.newElementSelected.bind(this)} />
                        <SelectedElements elements={this.state.selectedElements} onElementClicked={this.elementClicked.bind(this)} onElementDelete={this.deleteElement.bind(this)}/>
                    </div>

                    <div>
                        {alignResults}
                    </div>

                </CardText>

                <CardActions>
                    <FlatButton label="Query specified Elements" onClick={() => this.prepareResults()}/>
                </CardActions>
            </Card>

        );
    }
};

export interface ExplorePageProps { };
export interface ExplorePageState { queriesStored: number};

export class ExploreMainPage extends React.Component<ExplorePageProps, ExplorePageState> {

    allQueries = [];
    /**
     * Class constructor.
     */
    constructor(props) {
        super(props);

    }

    /**
     * This method will be executed after initial rendering.
     */
    componentWillMount() {

        this.setState({queriesStored: 0});

    }

    newQuery()
    {


        this.allQueries.push(<QueryComponent key={this.allQueries.length}/>);
        this.setState({queriesStored: this.allQueries.length});
        console.log("Query added");
    }

    clearQueries()
    {
        this.allQueries = [];
        this.setState({queriesStored: this.allQueries.length});
        console.log("Queries cleared");
    }

    /**
     * Render the component.
     */
    render() {

        console.log("ExploreMainpage render");


        return (

            <Paper>

                <Card>
                    <CardHeader
                    title="Without Avatar"
                    subtitle="Subtitle"
                    actAsExpander={true}
                    showExpandableButton={true}
                    />
                    <CardActions>
                    <FlatButton label="New Query" onClick={this.newQuery.bind(this)}/>
                    <FlatButton label="Clear" onClick={this.clearQueries.bind(this)}/>
                    </CardActions>
                    <CardText expandable={true}>

                        <p>Some explanation on how to use this view!</p>

                    </CardText>
                </Card>
                
                <div>
                {this.allQueries}
                </div>

            </Paper>

        );
    }

}