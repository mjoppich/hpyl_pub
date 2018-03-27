import * as React from 'react'; 
import {Card, CardActions, CardHeader, CardText} from 'material-ui/Card';
import FlatButton from 'material-ui/FlatButton';

import Paper from 'material-ui/Paper';

import SelectedElements from '../components/SelectedElements';
import ACInput from '../components/AutoComplete';
import MSATableViewer from '../components/MSATableViewer';

import axios from 'axios';
import config from '../config';

import OrganismChipAC from '../components/OrganismChipAC';
import EntityChipAC from '../components/EntityChipAC';

  
export interface QueryComponentProps { key: number};
export interface QueryComponentState { selectedElements: Array<any>, selectedOrganisms: Array<any>, alignments: any };

class QueryComponent extends React.Component<QueryComponentProps, QueryComponentState> {
    constructor(props) {
        super(props);

    }

    componentWillMount()
    {
        this.setState({selectedElements: []});
    }

    prepareResults()
    {
        var self = this;

        var allElems = [];
        for (var i = 0; i < this.state.selectedElements.length; ++i)
        {
            allElems.push(this.state.selectedElements[i].name);
        }

        axios.post(config.getRestAddress() + "/clustalign", {genes: allElems}, config.axiosConfig)
        .then(function (response) {
          console.log(response.data)

          self.setState({alignments: response.data})

        })
        .catch(function (error) {
          console.log(error)
          self.setState({alignments: {}})
        });


    }

    setSelectedEntities(elems : Array<any>)
    {
        console.log("Selected orgs:")
        console.log(elems)

        this.setState({selectedElements: elems});
    }

    setSelectedOrganisms(orgs : Array<any>)
    {
        console.log("Selected orgs:")
        console.log(orgs)

        this.setState({selectedOrganisms: orgs});
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

                    console.log("HOMALIGN");
                    console.log(homCluster);
                    alignResults.push( <MSATableViewer key={alignResults.length} alignments={homCluster}/> )
                }

            }
        }

        return (<Card style={{marginBottom: "20px"}}>
                <CardHeader
                title="Search Homology Entries"
                subtitle="Search by Gene/Protein ID"
                />
                <CardText>

                    <div>
                        
                        <EntityChipAC onValueChange={(entities) => this.setSelectedEntities(entities)}/>
                        <OrganismChipAC onValueChange={(allowedOrgs) => this.setSelectedOrganisms(allowedOrgs)}/>

                        <FlatButton label="Query specified Elements" onClick={() => this.prepareResults()}/>
                    </div>

                    <div>
                        {alignResults}
                    </div>

                </CardText>
            </Card>);

    //<ACInput onElementSelected={this.newElementSelected.bind(this)} />
    //<SelectedElements elements={this.state.selectedElements} onElementClicked={this.elementClicked.bind(this)} onElementDelete={this.deleteElement.bind(this)}/>
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

            <div>

                <Card style={{marginBottom: "20px"}}>
                    <CardHeader
                    title="Create Query"
                    subtitle="Manage your queries"
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

            </div>

        );
    }

}